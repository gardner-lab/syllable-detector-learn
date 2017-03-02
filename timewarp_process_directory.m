function [] = timewarp_process_directory(threshold, minsong, minnonsong, template_filename);
% Apply Nathan's dynamic timewarp to all the wave files in a directory.  I assume they're called
% "channel*.wav".  If the template file '0_template.wav' exists, it will be used, otherwise you will
% be prompted to create it.
%
% threshold:           match threshold below which stuff is considered a matching song (default all files)
% minsong, minnonsong: keep going until we exceed both of these numbers (default all files)
% template_filename:   filename of a wave file containing one sample of the target song (default 0_template.wav)
%
% OUTPUT:              creates a training file, song.mat, as specified in README.md
%
% DEPENDENCIES:        Nathan Perkins's find_audio, https://github.com/gardner-lab/find-audio

THRESHOLD_GAP = 1.5; % Nonsong must be above this factor of threshold.

files = dir('channel*.wav');
[~, sorted_index] = sortrows({files.date}');
files = files(sorted_index);
nfiles = length(files);

if ~exist('template_filename', 'var')
    template_filename = '0_template.wav';
end 
if ~exist('minsong', 'var')
    minsong = NaN;
end
if ~exist('minnonsong', 'var')
    minnonsong = NaN;
end


if isempty(template_filename) ...
        | ( ~exist(template_filename, 'file') & ~exist(strcat(template_filename, '.wav'), 'file'))
    make_template;
end
template = audioread(template_filename);
template_length = length(template);

song_n = 100;
song_i = 0;
song = zeros(template_length, song_n);

nonsong_n = 100;
nonsong_i = 0;
nonsong = zeros(template_length, nonsong_n);

if ~exist('wbar', 'var') | isempty(wbar) | ~ishandle(wbar) | isvalid(wbar)
    wbar = waitbar(0, 'Processing...');
else
    waitbar(0, wbar, 'Processing...');
end

for f = 1:nfiles
    waitbar(f/nfiles, wbar);
    
    if song_i >= minsong & nonsong_i >= minnonsong
        % Since the expensive step is find_audio and we have to run that for both song and non-song, there is no point in breaking
        % early if e.g. song_i >> minsong but nonsong_i < minnonsong.
        disp(sprintf('Reached the requested amount of data. %d songs, %d non-songs, %d%% of files processed.', ...
            song_i, nonsong_i, round(f*100/nfiles)));
        break;
    end
    
    [d, fs] = audioread(files(f).name);
    times = (1:length(d))/fs;
    allsamples_i = [];
    [starts, ends, scores] = find_audio(d, template, fs, 'threshold_score', threshold*THRESHOLD_GAP);
    disp(sprintf('%s: %d suspicious regions; %d above threshold', files(f).name, length(starts), length(find(scores <= threshold))));
    
    for n = 1:length(starts)
        % Toss them if they're not good enough
        if scores(n) > threshold
            continue;
        end
        sample_i = find(times >= starts(n) & times <= ends(n));
        allsamples_i = [allsamples_i sample_i]; % Keep track of these for later...
        sample = d(sample_i);
        sample_length = length(sample_i);
        sample_r = resample(sample, template_length, sample_length);
        
        % Preallocate some more memory, as required
        song_i = song_i + 1;
        if song_i > song_n
            song_n = song_n * 2;
            song(1, song_n) = 0;
        end
        
        song(:, song_i) = sample_r;
    end
    
    non_idx = setdiff(1:length(d), allsamples_i);
    nonsongdata = d(setdiff(1:length(d), allsamples_i));
    non_starter = 1;
    while non_starter + template_length - 1 < length(non_idx)
        nonsong_i = nonsong_i + 1;
        if nonsong_i > nonsong_n
            nonsong_n = nonsong_n * 2;
            nonsong(1, nonsong_n) = 0;
        end
        
        % Wrap the nonsong data into the appropriately sized chunks:
        nonsong(:, nonsong_i) = nonsongdata(non_starter:non_starter+template_length-1);
        non_starter = non_starter + template_length;
    end  
end
close(wbar);
wbar = [];

song = song(:, 1:song_i);
nonsong = nonsong(:, 1:nonsong_i);

save('song.mat', 'song', 'nonsong', 'fs');
