% Load data from an aggregate file (song.mat). Loads all the song data, and then only as much nonsong as is needed for the
% specified nonsinging_fraction.

% Copyright (C) 2017  Ben Pearre
%
% This file is part of the Zebra Finch Syllable Detector, syllable-detector-learn.
% 
% The Zebra Finch Syllable Detector is free software: you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at your option)
% any later version.
% 
% The Zebra Finch Syllable Detector is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with the Zebra Finch Syllable Detector.  If not, see
% <http://www.gnu.org/licenses/>.
function [ mic_data, spectrograms, nsamples_per_song, nmatchingsongs, nsongsandnonsongs, nonsinging_fraction, timestamps, nfreqs, freqs, ntimes, times, fft_time_shift_seconds, spectrogram_avg_img, power_img, freq_range_ds, time_window_steps, layer0sz, nwindows_per_song, noverlap] ...
    = load_roboaggregate_file(data_file, ...
    fft_time_shift_seconds_target, ...
    target_samplerate, ...
    fft_size, ...
    freq_range, ...
    time_window, ...
    nonsinging_fraction, ...
    n_whitenoise, ...
    song_crop_region)

load(data_file);

[orig_nsamples_per_song, nmatchingsongs] = size(song);

v = mean(var(song));

mic_data = song;

% Can we pull time-of-recording data out of the processing file?
timestamps = zeros(1, nmatchingsongs);
if exist('extract_filename', 'var')
    for i = 1:nmatchingsongs
        foo = strsplit(strrep(extract_filename{i}, '\', '/'), '/'); % Make sure directories are separated with /, and split
        cow = strsplit(foo{end}, '.'); % Filename component contains timestamp
        cor = strsplit(cow{6}, '-'); % Seems this field has a relevant number and then a -
        cow{6} = cor{1}; % reinsert expurgated value
        timestamps(i) = datenum([str2double(cow(2:6)) 0]); % pad seconds=0, store timestamp
    end
    
    disp(sprintf('File %s: %s -- %s', matching_song_file, datestr(timestamps(1)), datestr(timestamps(end))));
end




nnonmatches = size(nonsong, 2);

if nonsinging_fraction > 0
    if nnonmatches < nonsinging_fraction * nmatchingsongs
        nonsinging_fraction = nnonmatches/nmatchingsongs;
        warning('I had to lower nonsinging_fraction to %s.', sigfig(nonsinging_fraction));
    end
    if nnonmatches > nonsinging_fraction * nmatchingsongs
        nonsong = nonsong(:, 1:floor(nonsinging_fraction * nmatchingsongs));
    end
        
    mic_data = [mic_data nonsong];
end

    

%% Downsample the data to match target samplerate?
if fs ~= target_samplerate
        disp(sprintf('Resampling data from %g Hz to %g Hz...', fs, target_samplerate));
        [a b] = rat(target_samplerate/fs);
        
        mic_data = double(mic_data);
        mic_data = resample(mic_data, a, b);
end

%mic_data = mic_data / max(max(max(mic_data)), -min(min(mic_data)));

[nsamples_per_song, ~] = size(mic_data);

% Add some white noise.  This isn't really the right way to do this, but it may help:
mic_data = [mic_data wgn(nsamples_per_song, n_whitenoise, v, 'linear')];

[nsamples_per_song, nsongsandnonsongs] = size(mic_data);
song_seconds = nsamples_per_song / target_samplerate;
sample_times = 0:1/target_samplerate:song_seconds;
sample_times = sample_times(1:nsamples_per_song);

if ~isempty(song_crop_region)
    song_crop_i = find(sample_times >= song_crop_region(1)/1000 & sample_times <= song_crop_region(2)/1000);
    mic_data = mic_data(song_crop_i,:);
    % Repeat the above...
    [nsamples_per_song, nsongsandnonsongs] = size(mic_data);
end


%disp('Bandpass-filtering the data...');
%[B A] = butter(4, [0.03 0.9]);
%mic_data = single(filtfilt(B, A, double(mic_data)));


% Compute the spectrogram using original parameters (probably far from
% optimal but I have not played with them).  Compute one to get size, then
% preallocate memory and compute the rest in parallel.

noverlap = fft_size - (floor(target_samplerate * fft_time_shift_seconds_target));
% SPECGRAM(A,NFFT=512,Fs=[],WINDOW=[],noverlap=500)
%speck = specgram(mic_data(:,1), 512, [], [], 500) + eps;

window = hamming(fft_size);

[speck freqs times] = spectrogram(mic_data(:,1), window, noverlap, [], target_samplerate);
% Adjust "times" to reflect the time at which the information is actually available--i.e. the end,
% rather than the middle, of the window:
times = times - times(1) + fft_size/target_samplerate;

[nfreqs, ntimes] = size(speck);
speck = speck + eps;

% This will be approximately the same as fft_time_shift_seconds_target, but not quite: the fft_time_shift
% is given by noverlap, and will actually be fft_size/target_samplerate
fft_time_shift_seconds = (times(end)-times(1))/(length(times)-1);
fprintf('FFT time shift = %d frames, %s... ms\n', ...
    floor(target_samplerate * fft_time_shift_seconds_target), ...
    sigfig(1000*fft_time_shift_seconds, 8));


spectrograms = zeros([nsongsandnonsongs nfreqs ntimes]);
spectrograms(1, :, :) = speck;
disp('Computing spectrograms...');
parfor i = 2:nsongsandnonsongs
        spectrograms(i, :, :) = spectrogram(mic_data(:,i), window, noverlap, [], target_samplerate) + eps;
end

spectrograms = single(spectrograms);


% Create a pretty graphic for display (which happens later)
spectrograms = abs(spectrograms);
spectrogram_avg_img = squeeze(sum(spectrograms(1:nmatchingsongs,:,:)));

% Poor man's filter to take off the high amplitudes at uselessly low frequencies:
low_freq_i = find(freqs < 200);
spectrogram_avg_img(low_freq_i,:) = 0;

% Kill off any signal below mean-1std, making the spectrogram prettier:
spectrogram_avg_img = reshape(zscore(spectrogram_avg_img(:)), size(spectrogram_avg_img)) + 1;
spectrogram_avg_img(find(spectrogram_avg_img < 0)) = 0;
spectrogram_avg_img = log(spectrogram_avg_img + eps);

%% Draw the pretty full-res spectrogram
try
    figure(4);
    imagesc([times(1) times(end)]*1000, [freqs(1) freqs(end)]/1000, spectrogram_avg_img);
    axis xy;
    xlabel('Time (ms)');
    ylabel('Frequency (kHz)');
end

% Construct "ds" (downsampled) dataset.  This is heavily downsampled to save on computational
% resources.  This would better be done by modifying the spectrogram's
% parameters above (which would only reduce the number of frequency bins,
% not the number of timesteps), but this will do for now.

% Number of samples: (nsongs*(ntimes-time_window))
% Size of each sample: (ntimes-time_window)*length(freq_range)

%%%%%%%%%%%%


freq_range_ds = find(freqs >= freq_range(1) & freqs <= freq_range(2));
disp(sprintf('Using frequencies in [ %g %g ] Hz: %d frequency samples.', ...
    freq_range(1), freq_range(2), length(freq_range_ds)));
time_window_steps = double(floor(time_window / fft_time_shift_seconds));
disp(sprintf('Time window is %g ms, %d samples.', time_window*1000, time_window_steps));

% How big will the neural network's input layer be?
layer0sz = length(freq_range_ds) * time_window_steps;

% The training input set X is made by taking all possible time
% windows.  How many are there?  The training roboaggregateput set Y will be made by
% setting all time windows but the desired one to 0.
nwindows_per_song = ntimes - time_window_steps + 1;

disp('Creating spectral power image...');

% Create an image on which to superimpose the results...
power_img = squeeze((sum(spectrograms, 2)));
power_img(find(isinf(power_img))) = 0;
