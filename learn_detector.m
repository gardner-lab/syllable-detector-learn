% learn_detector: train a neural network to detect zebra finch syllables.
%
%        Requires data in (by default) 'song.mat', and training configuration
%        in 'params.m' and/or parameters given as 
%        learn_detector('parameter', value) pairs. See README.md for
%        instructions.
%
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
    
function learn_detector(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% These are defaults.  If you want to change them, do so in params.m
times_of_interest_ms = NaN;                      % Milliseconds at which to trigger. Exists here just to define it as a valid parameter.
nhidden_per_output = 4;                          % How many hidden units per syllable?  2 works and trains fast.  4 works ~20% better...
fft_time_shift_seconds_target = 0.0015;          % FFT frame rate (seconds).  Paper mostly used 0.0015 s: great for timing, but slow to train.
use_jeff_realignment_train = false;              % Micro-realign at each detection point using Jeff's time-domain code?  Don't do this.
use_jeff_realignment_test = false;               % Micro-realign test data only at each detection point using Jeff's time-domain code.  Nah.
use_nn_realignment_test = false;                 % Try using the trained network to realign test songs (reduce jitter?)
nonsinging_fraction = 1;                         % Train on this proportion of nonsinging data (e.g. cage noise, calls)
n_whitenoise = 10;                               % Add this many white noise samples (may help with white noise stimulation, but not as much as adding real samples!
ntrain_approx_max_matching_songs = 1000;         % Total dataset size is songs+nonsongs.  Only this+this*nonsinging_fraction will be used to train, leaving the rest for test
testfile_include_nonsinging = false;             % Include nonsinging data in audio test file (no point if just used to measure timing)
samplerate = 44100;                              % Target samplerate should match sampling frequency of live detector
fft_size = 256;                                  % FFT size
use_pattern_net = false;                         % Use MATLAB's pattern net (fine, but no control over false-pos vs false-neg cost)
do_not_randomise = false;                        % Use songs in original order?
separate_network_for_each_syllable = true;       % Train a separate network for each time of interest?  Or one network with multiple outs?
nruns = 1;                                       % Perform a few training runs and create beeswarm plot? (paper figure 3 used 100)
freq_range = [1000 8000];                        % Frequencies of the song to examine.
time_window_ms = 50;                             % How many seconds long is the detection sliding-time-window? Useful range is about 30-100.
false_positive_cost = 1;                         % Cost of false positives is relative to that of false negatives.
create_song_test_file = -1;                      % Big, and take some time to save.  1: always, 0: never, -1: only if nruns == 1
song_crop_region = [-Inf Inf];                   % Crop the ends off the aligned song (milliseconds).  [] or [-Inf Inf] for no crop.  YOU MAY (OR MAY NOT?) WANT TO ADJUST TIMES OF INTEREST: times_of_interest_ms = [ 550 660 ] - song_crop_region(1)
log_file_exists_action = 'ask';                  % If confusion_log exists, 'ask', 'append', 'replace'


% The two required files:
data_file = 'song';                              % song.mat
params_file = 'params';                          % params.m



%use_previously_trained_network = '5syll_1ms.mat' % Rather than train a new network, use this one? NO ERROR CHECKING!!!!!
%  Finally: where do the aligned song and nonsong data files live?  And which times do we care
%  about?

% Load the user configuration.  This is done by a function that runs the params file as a .m file and adds all the discovered
% parameters to the struct 'p'. This can then be checked against variables that actually exist.
if ~isempty(params_file)
    disp(sprintf('********** Configuration file %s: *************', ...
        strcat(pwd, filesep, params_file, '.m')));
    
    user_parameters = load_params(params_file)
    
    pf = fieldnames(user_parameters);
    for i = 1:length(pf)
        if exist(pf{i}, 'var')
            eval(sprintf('%s = user_parameters.%s;', pf{i}, pf{i}));
        else
            error('Parameter name ''%s'' is invalid.', pf{i});
        end
    end
end

if nargin > 0
    disp('********** Processing function arguments... ******************************');
    
    % Also allow override of any parameter on the command line
    for i = 1:2:nargin
        if ~exist(varargin{i}, 'var')
            warning('learn_detector: Argument ''%s'' is invalid.', varargin{i});
        else
            sprintf('... %s', varargin{i});
            eval(sprintf('%s = varargin{i+1};', varargin{i}));
        end
    end
    
    disp('**************************************************************************');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% End Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltasong = false;
if deltasong % Delta function is a special test case for measuring detector latency.
    agg_audio.fs = 44100;
    indices = round(-0.010 * agg_audio.fs);
    times_of_interest_ms = 300;
    samples_of_interest = round((times_of_interest_ms / 1e3) * agg_audio.fs) + 1;
    n = 128;
    mic_data = rand([20000, n])/100;
    mic_data(samples_of_interest + indices, :) = rand([length(indices), n])/100 + 1;
end


if exist('use_previously_trained_network', 'var') & ~isempty(use_previously_trained_network)
    disp(sprintf('Loading previously trained network ''%s''...', use_previously_trained_network));
    load(use_previously_trained_network);
end


rng('shuffle');



% If confusion_log_perf.txt exists, there is the risk that something important (parameters, code...)
% has changed since that file was last added to.  Ask the user.
if nruns > 1 & separate_network_for_each_syllable & exist('./confusion_log_perf.txt', 'file')
    if strcmp(log_file_exists_action, 'ask')
        log_file_exists_action = questdlg('A multi-run logfile exists.  Keep adding to it?', ...
            'What should I do with the logfile?', ...
            'append', ...
            'replace', ...
            'append');
    end
    
    if strcmp(log_file_exists_action, 'replace')
        old_file = dir('confusion_log_perf.txt');
        old_file_date = datestr(old_file.datenum, 31);
        old_file_new_name = sprintf('confusion_log_perf.%s.txt', old_file_date);
        
        movefile('confusion_log_perf.txt', old_file_new_name, 'f');
        first_run = 1;
        disp(sprintf('Old logfile has been renamed ''%s''', old_file_new_name));
    elseif strcmp(log_file_exists_action, 'append')
        % This file is created in show_confusion.m.
        confusion = load('confusion_log_perf.txt');
        %sylly = union(confusion(:,1)', times_of_interest_ms / 1000));
        % But we probably only care about the ones in times_of_interest_ms. So scratch that?
        sylly = times_of_interest_ms/1000;
        sylly_counts = [];
        for i = 1:length(sylly)
            sylly_counts(i) = length(find(confusion(:,1)==sylly(i)));
        end
        syllables_still_needed = [sylly * 1000 ; nruns - sylly_counts];
        if ~sum(syllables_still_needed(2,:))
            disp('Done! (Looks like we have all the trials we needed.)');
            return;
        end

        disp(sprintf('Continuing where we left off. Progress so far (time_of_interest_ms ; # runs still required):'));
        syllables_still_needed
        
        first_run = min(sylly_counts) + 1;
    else
        error('log_file_exists_action = ''%s'', but must be ''append'', ''replace'', or ''ask''', log_file_exists_action);
    end
else
    first_run = 1;
end

time_window_s = time_window_ms / 1e3;


[ mic_data, spectrograms, nsamples_per_song, nmatchingsongs, nsongsandnonsongs, nonsinging_fraction, timestamps, ...
    nfreqs, freqs, ntimes, times, fft_time_shift_seconds, spectrogram_avg_img_songs_log, spectrogram_power_img, ...
    freq_range_ds, time_window_steps, layer0sz, nwindows_per_song, noverlap] ...
    = load_roboaggregate_file(data_file, ...
    fft_time_shift_seconds_target, ...
    samplerate, ...
    fft_size, ...
    freq_range, ...
    time_window_s, ...
    nonsinging_fraction, ...
    n_whitenoise, ...
    song_crop_region);

fprintf('Got %d songs, %d cage noise, %d songs-and-nonsongs including %d of synthetic white noise.\nNonsinging fraction is %s.\n', ...
    nmatchingsongs, nsongsandnonsongs-nmatchingsongs-n_whitenoise, nsongsandnonsongs, n_whitenoise, sigfig(nonsinging_fraction, 4));


%% Draw the spectral image.  If no times_of_interest defined, this is what the user will use to choose some.
try
    figure(4);
    subplot(1,1,1);
    specfig = imagesc(times([1 end])*1000, freqs([1 end])/1000, spectrogram_avg_img_songs_log);
    axis xy;
    xlabel('Time (ms)');
    ylabel('Frequency (kHz)');
    set(gca, 'YLim', [0 10]);
end
if ~exist('times_of_interest_ms', 'var') | isempty(times_of_interest_ms)
    disp(sprintf('No times of interest defined.  Please look at the spectrogram in Figure 4 and define one or more in ''%s'', with "times_of_interest_ms = [x y];" for detection at x and y milliseconds into the spectrogram.', strcat(pwd, filesep, params_file, '.m')));
    return;
end

times_of_interest_s = times_of_interest_ms / 1e3;



%% Define training set
% Hold some data out for final testing.  This includes both matching and non-matching IF THE SONGS
% ARE IN RANDOM ORDER
ntrainsongsandnonsongs = min(floor(nsongsandnonsongs*8/10), floor((1+nonsinging_fraction)*ntrain_approx_max_matching_songs));
ntestsongsandnonsongs = nsongsandnonsongs - ntrainsongsandnonsongs;
disp(sprintf('%d training songs-and-nonsongs.  %d remain for test.', ntrainsongsandnonsongs, ntestsongsandnonsongs));

% If we're using "fit", it'll produce useless warnings (some kludgey analysis I do later uses "fit",
% but I want to disable them outside the loop).  Silence them!
warning('off', 'curvefit:prepareFittingData:nonDouble');
warning('off', 'curvefit:prepareFittingData:sizeMismatch');
warning('off', 'curvefit:prepareFittingData:removingNaNAndInf')


% Just one rudimentary error-check:
if any(times_of_interest_ms < time_window_ms) | any(times_of_interest_ms > times(end) * 1000)
    error('learn_detector:invalid_time', ...
        'All times_of_interest_ms [ %s] must be >= time_window_ms (%g) and < %s', ...
        sprintf('%g ', times_of_interest_ms), time_window_ms, times(end) * 1000);
end

% Create informative names for the detection points:
if ~exist('times_of_interest_names', 'var') | length(times_of_interest_names) < length(times_of_interest_separate)
    for i = 1:length(times_of_interest_ms)
        times_of_interest_names{i} = sprintf('t^*_{%d}', round(times_of_interest_ms(i)));
    end
end



% Create a FOR loop over these, if necessary
if separate_network_for_each_syllable
    times_of_interest_separate = times_of_interest_s;
else
    times_of_interest_separate = NaN;
    times_of_interest_simultaneous = times_of_interest_s;
end

training_times = [];
loop_times = [];
catch_up = false;

tic;
start_time = datetime('now');
eta = 'next weekend';


for run = first_run:nruns
    disp(sprintf('Starting run #%d...', run));
    
    
        
    separate_syllable_counter = 0;
    
    for thetime = times_of_interest_separate
        
        % When we've completed only a proper subset of times_of_interest within a run, do the run but skip the redundant
        % times_of_interest.
        if separate_network_for_each_syllable ...
                & exist('sylly_counts', 'var')
            n_completed_of_this_time_of_interest = sylly_counts(find(sylly == thetime));
            if run <= n_completed_of_this_time_of_interest 
                disp(sprintf('We already finished timepoint %g for run %d. Continuing...', thetime, run));
                continue;
            end
        end
        % thetime will be each of the times_of_interest, or else NaN, which will run once through the
        % loop.
        
        separate_syllable_counter = separate_syllable_counter + 1;


        % On each run of this loop, change the presentation order of the
        % data, so we get (a) a different subset of the data than last time for
        % training vs. final testing and (b) different training data presentation
        % order.
        rng('shuffle');
        
        if do_not_randomise
            randomorder(1:nsongs) = 1:nsongsandnonsongs;
            warning('NOT permuting song order');
        else
            % Because load_roboaggregate_file() makes sure that nonsinging_fraction is accurate, drawing randomly from
            % songs-and-nonsongs will give roughly the correct ratio of song to nonsong.
            randomorder = randperm(nsongsandnonsongs);
        end
                
        if separate_network_for_each_syllable
            % "toi" will be times_of_interest_s(separate_syllable_counter)
            toi = thetime;
        else
            % This is redundant, but here for readability:
            toi = times_of_interest_s;
        end
        
        
        
        ntsteps_of_interest = length(toi);
        for i = 1:length(toi)
            tsteps_of_interest(i) = find(times >= toi(i), 1);
        end
        
        disp(sprintf('********** Working on [ %s] ms **********', sprintf('%g ', toi * 1000)));
        
        
        %% Create the training set
        if deltasong
            shotgun_sigma = 0.00001;
        else
            shotgun_sigma = 0.002; % TUNE
        end
        
        
        
        
        
        if use_jeff_realignment_train | use_jeff_realignment_test
            
            %% For each timestep of interest, get the offset of this song from the most typical one.
            disp('Computing target jitter compensation...');
            
            % We'll look for this long around the timestep, to compute the canonical
            % song
            time_buffer = 0.04;
            tstep_buffer = round(time_buffer / fft_time_shift_seconds);
            
            % For alignment: which is the most stereotypical song at each target?
            
            %[B A] = butter(4, [0.01 0.05]);
            %mic_data2 = filtfilt(B, A, double(mic_data));
            
            for i = 1:ntsteps_of_interest
                range = tsteps_of_interest(i)-tstep_buffer:tsteps_of_interest(i)+tstep_buffer;
                range = range(find(range>0&range<=ntimes));
                foo = reshape(spectrograms(1:nmatchingsongs, :, range), nmatchingsongs, []) * reshape(mean(spectrograms(:, :, range), 1), 1, [])';
                [val canonical_songs(i)] = max(foo);
                [target_offsets(i,:) sample_offsets(i,:)] = get_target_offsets_jeff(mic_data(:, 1:nmatchingsongs), tsteps_of_interest(i), samplerate, fft_time_shift_seconds, canonical_songs(i));
            end
            
            target_offsets_test = target_offsets;
            sample_offsets_test = sample_offsets;
            if ~use_jeff_realignment_train
                fprintf('\n               ***** DISCARDING TARGET JITTER COMPENSATION FOR TRAINING *****\n\n');
                target_offsets = 0 * target_offsets;
                sample_offsets = 0 * sample_offsets;
            end
            if ~use_jeff_realignment_test
                fprintf('\n               ***** DISCARDING TARGET JITTER COMPENSATION FOR TEST FILE *****\n\n');
                target_offsets_test = 0 * target_offsets_test;
                sample_offsets_test = 0 * sample_offsets_test;
            end
        else
            target_offsets = zeros(ntsteps_of_interest, nsongsandnonsongs);
            sample_offsets = target_offsets;
            target_offsets_test = target_offsets;
            sample_offsets_test = sample_offsets;
        end
        %hist(target_offsets', 40);
        
        %pn = 1:nmatchingsongs;
        %[vt pt] = sort(target_offsets);
        %[vs ps] = sort(sample_offsets);
        %figure(4);
        %subplot(1,1,1);
        %power_img = power_img(1:nmatchingsongs,:);
        %imagesc(power_img(pt,:));
        %set(gca, 'xlim', [280.2 300]);
        
        
        
        %% Draw the pretty full-res spectrogram and the targets
        try
            figure(4);
            subplot(1,1,1);
            %subplot(ntsteps_of_interest+1,1,1);
            specfig = imagesc(times([1 end])*1000, freqs([1 end])/1000, spectrogram_avg_img_songs_log);
            axis xy;
            xlabel('Time (ms)');
            ylabel('Frequency (kHz)');
            
            % Draw the syllables of interest:
            for i = 1:ntsteps_of_interest
                line(toi(i)*[1;1]*1000, freqs([1 end])/1000, 'Color', [1 0 0]);
                windowrect = rectangle('Position', [(toi(i) - time_window_s)*1000 ...
                    freq_range(1)/1000 ...
                    time_window_s(1)*1000 ...
                    (freq_range(2)-freq_range(1))/1000], ...
                    'EdgeColor', [1 0 0]);
            end
            set(gca, 'YLim', [0 10]);
            drawnow;
        catch ME
        end
        
        disp(sprintf('Creating training set from %d songs and nonsongs; test set from the remainder...', ntrainsongsandnonsongs));
        % These are shuffled according to randomorder. Because nonsingin_fraction is used in the creation of the data set, any
        % random set of these (e.g. 1:1000) will have the requested mix of song:nonsong.
        [nnsetX nnsetY] = create_training_set(spectrograms, ...
            tsteps_of_interest, ...
            target_offsets, ...
            shotgun_sigma, ...
            randomorder, ...
            nmatchingsongs, ...
            nsongsandnonsongs, ...
            nwindows_per_song, ...
            layer0sz, ...
            fft_time_shift_seconds, ...
            time_window_steps, ...
            ntimes, ...
            freq_range_ds);
        
        if use_pattern_net
            nnsetYC = [nnsetY~=0 ; nnsetY==0];
        end
        
        %yy=reshape(nnsetY, nwindows_per_song, nsongsandnonsongs);
        %imagesc(yy');
        
        % original order: spectrograms, spectrograms_ds, song_montage
        %   indices into original order: trainsongs, testsongs
        % shuffled: nnsetX, nnsetY, testout
        %   indices into shuffled arrays: nnset_train, nnset_test
        
        % We can make these contiguous blocks, since the spectrograms were
        % shuffled when the training sets were built.
        nnset_train = 1:(ntrainsongsandnonsongs * nwindows_per_song);
        nnset_test = ntrainsongsandnonsongs * nwindows_per_song + 1 : size(nnsetX, 2);
        
        % Create the network.  The parameter is the number of units in each hidden
        % layer.  [8] means one hidden layer with 8 units.  [] means a simple
        % perceptron.
        
        if use_pattern_net
            net = patternnet(nhidden_per_output * ntsteps_of_interest);
        else
            net = feedforwardnet([nhidden_per_output * ntsteps_of_interest]);
        end
        net.inputs{1}.processFcns={'mapstd'};
        %net.trainFcn = 'trainlm';
        net.trainFcn = 'trainscg';
        
        
        net.divideParam.trainRatio = 80/100;
        net.divideParam.valRatio = 20/100;
        net.divideParam.testRatio = 0/100;
        net.plotFcns = {'plotperform'};
        %net.trainParam.goal=1e-3;
        
        fprintf('(%s) Training network with %s...\n', datestr(datetime('now'), 'HH:MM'), net.trainFcn);
        
        
        % Once the validation set performance stops improving, it seldom seems to
        % get better, so keep this small.  OOPS--that was true for trainlm, but trainscg has much
        % faster (and more) iterations.
        %net.trainParam.max_fail = 3;
        
        loop_times(end+1) = toc/60;
        
        tic
        if exist('use_previously_trained_network', 'var') & ~isempty(use_previously_trained_network)
            load(use_previously_trained_network);
        else
            if use_pattern_net
                [net, train_record] = train(net, nnsetX(:, nnset_train), nnsetYC(:, nnset_train));
            else
                [net, train_record] = train(net, nnsetX(:, nnset_train), nnsetY(:, nnset_train));
            end
        end
        
        % Oh yeah, the line above was the hard part.
        training_times = [training_times toc/60];
        if length(loop_times) > 1
            disp(sprintf('   ...training took %g minutes (mean %s m).  Rest-of-loop %g (mean %s m).', ...
                toc/60, sigfig(mean(training_times)), ...
                loop_times(end), sigfig(mean(loop_times(2:end)))));
            %length(times_of_interest_ms) * (nruns-run) * mean(loop_times(2:end))
        end
        
        
        tic
        % Test on all the data:
        
        % Test just on the non-training data, right?  Compute them all, and then only count ntestsongsandnonsongs for statistics (later)
        if use_pattern_net
            testout = net(nnsetX);
            testout = testout(1,:);
        else
            testout = sim(net, nnsetX);
        end
        testout = reshape(testout, ntsteps_of_interest, nwindows_per_song, nsongsandnonsongs);
        
        % Update the each-song image
        power_img = spectrogram_power_img(randomorder,:);
        power_img = repmat(power_img / max(max(power_img)), [1 1 3]);
        
        disp('Computing optimal output thresholds...');
        
        % How many seconds on either side of the tstep_of_interest is an acceptable match?
        MATCH_PLUSMINUS = 0.02;
        
        % Which songs should have hits?  The first nmatchingsongs, but permuted to the same order as the
        % training/test sets, as given by randomorder.
        songs_with_hits = [ones(1, nmatchingsongs) zeros(1, nsongsandnonsongs - nmatchingsongs)]';
        songs_with_hits = songs_with_hits(randomorder);
        
        % Search for the optimal trigger thresholds using just the training set
        trigger_thresholds = optimise_network_output_unit_trigger_thresholds(...
            testout(:,:,1:ntrainsongsandnonsongs), ...
            nwindows_per_song, ...
            false_positive_cost, ...
            toi, ...
            tsteps_of_interest, ...
            MATCH_PLUSMINUS, ...
            fft_time_shift_seconds, ...
            time_window_steps, ...
            songs_with_hits(1:ntrainsongsandnonsongs), ...
            true);
        
        % Now that we've computed the thresholds using just the training set, print the accuracy (as confusion matrices)
        % using just the holdout test set.
        foo = ntrainsongsandnonsongs+1:size(testout, 3);
        show_confusion(...
            testout(:, :, foo), ...
            nwindows_per_song, ...
            false_positive_cost, ...
            toi, ...
            tsteps_of_interest, ...
            MATCH_PLUSMINUS, ...
            fft_time_shift_seconds, ...
            time_window_steps, ...
            songs_with_hits(foo), ...
            trigger_thresholds);
        
        
        %figure(32);
        %plot(times(time_window_steps:end), squeeze(testout(1,:,:)), 'b', ...
        %    times([time_window_steps end]), [1 1]*trigger_thresholds, 'r');
        %title('Network output and threshold');
        
        SHOW_THRESHOLDS = true;
        SHOW_ONLY_TRUE_HITS = true;
        SORT_BY_ALIGNMENT = true;
        raster_colour_left_bar = false;
        % For each fft_time_shift_seconds of interest, draw that output unit's response to all
        % timesteps for all songs:
                
        target_offsets_net = zeros(ntsteps_of_interest, nsongsandnonsongs);
        sample_offsets_net = zeros(ntsteps_of_interest, nsongsandnonsongs);

        % For each TOI, plot its response graph
        try
        figure(6);
        for i = 1:length(toi)
            if separate_network_for_each_syllable
                subplot(length(times_of_interest_separate), 1, separate_syllable_counter);
            else
                subplot(ntsteps_of_interest, 1, i);
            end
            
            testout_i_squeezed = reshape(testout(i,:,:), [], nsongsandnonsongs);
            leftbar = zeros(time_window_steps-1, nsongsandnonsongs);
            
            if SHOW_THRESHOLDS
                % "img" is a tricolour image
                img = power_img;
                % de-bounce:
                trigger_img = trigger(testout_i_squeezed', trigger_thresholds(i), 0.1, fft_time_shift_seconds);
                trigger_img = [leftbar' trigger_img];
                [val pos] = max(trigger_img, [], 2);
                triggertimes(:,i) = pos * fft_time_shift_seconds;
                
                [targets_with_offsets, target_offsets_net_tmp] = find(trigger_img);
                target_offsets_net(i,targets_with_offsets) = target_offsets_net_tmp' - tsteps_of_interest(i) + 1;
                sample_offsets_net(i,:) = target_offsets_net(i,:) * fft_time_shift_seconds * samplerate;
             
                %figure(7);
                %hist([target_offsets_2 ; target_offsets_net]', 50);
                %hist([sample_offsets_test ; sample_offsets_net]', 50);
                %target_offset_mean_difference = mean(target_offsets_test) - mean(target_offsets_net)
                %figure(6);
                
                % img is RGB.  Here I'm playing with colouring the image with triggers
                img(1:ntrainsongsandnonsongs, :, 1) = img(1:ntrainsongsandnonsongs, :, 1) - trigger_img(1:ntrainsongsandnonsongs, :);
                img(1:ntrainsongsandnonsongs, :, 2) = img(1:ntrainsongsandnonsongs, :, 2) + trigger_img(1:ntrainsongsandnonsongs, :);
                img(1:ntrainsongsandnonsongs, :, 3) = img(1:ntrainsongsandnonsongs, :, 3) + trigger_img(1:ntrainsongsandnonsongs, :);
                % Different colour for testsongs
                img(ntrainsongsandnonsongs+1:end, :, 1) = img(ntrainsongsandnonsongs+1:end, :, 1) + trigger_img(ntrainsongsandnonsongs+1:end, :);
                img(ntrainsongsandnonsongs+1:end, :, 2) = img(ntrainsongsandnonsongs+1:end, :, 2) - trigger_img(ntrainsongsandnonsongs+1:end, :);
                img(ntrainsongsandnonsongs+1:end, :, 3) = img(ntrainsongsandnonsongs+1:end, :, 3) - trigger_img(ntrainsongsandnonsongs+1:end, :);
                
                if raster_colour_left_bar
                    % Colour the leftbar according to train and test:
                    img(1:ntrainsongsandnonsongs, 1:time_window_steps, 3) = 1;
                    img(1:ntrainsongsandnonsongs, 1:time_window_steps, 2) = 1;
                    img(1:ntrainsongsandnonsongs, 1:time_window_steps, 1) = 0;
                    img(ntrainsongsandnonsongs+1:end, 1:time_window_steps, 2) = 0;
                    img(ntrainsongsandnonsongs+1:end, 1:time_window_steps, 1) = 1;
                    img(ntrainsongsandnonsongs+1:end, 1:time_window_steps, 3) = 0;
                end
                
                if SHOW_ONLY_TRUE_HITS
                    img = img(find(songs_with_hits), :, :);
                    pos = pos(find(songs_with_hits));
                end
                
                if SORT_BY_ALIGNMENT
                    %[a, new_world_order] = sort(sample_offsets(randomorder(1:nmatchingsongs)));
                    [~, new_world_order] = sort(pos);
                    img = img(new_world_order,:,:);
                end
                
                % Make sure the image handle has the correct axes
                if SHOW_ONLY_TRUE_HITS
                    imh = image([times(1) times(end)]*1000, [1 sum(songs_with_hits)], img);
                else
                    imh = image([times(1) times(end)]*1000, [1 nsongsandnonsongs], img);
                end
            else
                leftbar(:, 1:ntrainsongsandnonsongs) = max(max(testout_i_squeezed))/2;
                leftbar(:, ntrainsongsandnonsongs+1:end) = 3*max(max(testout_i_squeezed))/4;
                testout_i_squeezed = [leftbar' testout_i_squeezed'];
                imagesc([times(1) times(end)]*1000, [1 nsongsandnonsongs], testout_i_squeezed);
            end
            xlabel('Time (ms)');
            ylabel('Song');
            %title(sprintf('Detection events for %d ms', round(1000*toi(i))));
            if ~catch_up
                title(sprintf('Detection events: %s', times_of_interest_names{separate_syllable_counter}));
            end
            
            if ~SORT_BY_ALIGNMENT
                %% Show coloration by labeling the blocks of training and test songs
                text(time_window_s/2*1000, ntrainsongsandnonsongs/2, 'train', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
                text(time_window_s/2*1000, ntrainsongsandnonsongs+ntestsongsandnonsongs/2, 'test', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
            elseif raster_colour_left_bar
                %% Show colouration by labeling the largest contiguous blocks of training and test songs
                s = img(:,1,1); % s is now 0 for train, 1 for test
                a = diff(s);
                b = find([a; Inf] ~= 0);
                c = diff([0; b]);
                d = cumsum(c);
                [e f] = max(c(2:2:end)); % even: test
                f = f * 2;
                testcentre = mean(d(f-1:f));
                
                [g h] = max(c(1:2:end)); % odd: train
                h = h * 2 - 1;
                try
                    traincentre = mean(d(h-1:h));
                catch ME
                    traincentre = mean([0 d(h)]);
                end
                
                if s(1) % Fix parity if it looks to be wrong...
                    foo = testcentre;
                    testcentre = traincentre;
                    traincentre = foo;
                end
                text(time_window_s/2*1000+3, traincentre, 'train', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 0);
                text(time_window_s/2*1000+3, testcentre, 'test', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 0);
            else
                %% Show colouration via a legend
                if separate_syllable_counter == 1
                    prevhold = ishold;
                    hold on;
                    % Plot some dummy lines offscreen for legend to pick up the colours from:
                    xlims = get(gca, 'XLim');
                    ylims = get(gca, 'YLim');
                    plot([-1 -2], [-1 -2], 'color', [0 1 1], 'LineWidth', 5);
                    plot([-1 -2], [-1 -2], 'color', [1 0 0], 'LineWidth', 5);
                    if ~prevhold
                        hold off;
                    end
                    set(gca, 'XLim', xlims, 'YLim', ylims);
                    
                    legend('Train', 'Test', 'location', 'SouthWest');
                end
            end
        end
        catch ME
        end
        
        %% If possible, plot variability over the course of the day.  This
        % requires one network trained on multiple syllables, and computes timing differences
        % between each of the syllables in the roboaggregate file.
        if ntsteps_of_interest >= 2 & false
            triggertimes_backup = triggertimes;
            
            % For syllable timing variability analysis, look only for "correct" hits.  Is this what I want?
            triggertimes(find(triggertimes == 0)) = NaN;
            for i = 1:ntsteps_of_interest
                ntt = triggertimes(:,i);
                ntt(find(abs(ntt - toi(i)) > 0.02)) = NaN;
                triggertimes(:,i) = ntt;
            end
            
            %figure(233);
            %clf;
            %hold on;
            
            tsn = (timestamps - timestamps(1))*24;           
            tsu = unique(tsn);
            
            ncombs = nchoosek(ntsteps_of_interest, 2);
            combs = nchoosek(1:ntsteps_of_interest, 2);
            colours = distinguishable_colors(ncombs);
            
            comb = 0;
            for first = 1:ntsteps_of_interest
                for second = first+1:ntsteps_of_interest
                    comb = comb + 1;
                    for i = 1:length(tsu)
                        inds = find(tsn(randomorder(1:nmatchingsongs)) == tsu(i));
                        measure = (triggertimes(inds, first) - triggertimes(inds, second)) * 1e3;
                        binsize(i,first,second) = length(measure) - sum(isnan(measure));
                        means(i,first,second) = nanmean(measure);
                        stds(i,first,second) = nanstd(measure);
                        %plot(tsu(i)*ones(size(measure)), measure, '.', 'Color', colours(comb,:));
                    end
                end
            end
            ste = stds ./ binsize.^(1/2);
            %hold off;
            stds(find(stds==0)) = NaN;
            
            figure(234);
            clf;
            for first = 1:ntsteps_of_interest
                for second = first+1:ntsteps_of_interest
                    subplot(ntsteps_of_interest-1, ntsteps_of_interest-1, (first-1)*(ntsteps_of_interest-1)+second-1);
                    includ = find(~isnan(stds(:,first,second)));
                    [tsup stdsp] = prepareCurveData(tsu, stds(:,first,second));
                    xl = [min(tsup) max(tsup)];
                    if length(tsup) > 2
                        f0 = fit(tsup, stdsp, 'poly1', 'weights', binsize(includ,first,second));
                        hold on;
                        scatter(tsup, stdsp, binsize(includ,first,second), [0 0 1], 'filled');
                        plot(xl, f0(xl), 'Color', [0 0 1]);
                        hold off;
                        xlabel('Time (hours)');
                        ylabel('Intersyllable std dev (ms)');
                        yl = get(gca, 'YLim');
                        yl(1) = 0;
                        set(gca, 'XLim', xl, 'YLim', yl);
                        confs = confint(f0);
                        title(sprintf('%d-%d: slope %s (95%%)', first, second, sigfig(confs(1:2), 2)));
                    end
                end
            end
        end
        
        drawnow;
        
        
        
        %% Realign test output according to the neural network's detection point?
        if use_nn_realignment_test
            target_offsets_net = zeros(ntsteps_of_interest, nsongsandnonsongs);
            sample_offsets_net = zeros(ntsteps_of_interest, nsongsandnonsongs);
            for i = 1:ntsteps_of_interest
                testout_i_squeezed = reshape(testout(i,:,:), [], nsongsandnonsongs);
                leftbar = zeros(time_window_steps-1, nsongsandnonsongs);
                trigger_img = trigger(testout_i_squeezed', trigger_thresholds(i), 0.1, fft_time_shift_seconds);
                trigger_img = [leftbar' trigger_img];
                [val pos] = max(trigger_img, [], 2);
                
                [targets_with_offsets, target_offsets_net_tmp] = find(trigger_img);
                target_offsets_net(i,targets_with_offsets) = target_offsets_net_tmp' - tsteps_of_interest(i) + 1;
                sample_offsets_net(i,:) = target_offsets_net(i,:) * fft_time_shift_seconds * samplerate;
                target_offsets_test = target_offsets_net;
                sample_offsets_test = sample_offsets_net;
            end
        end
        
        
        if nruns > 1 & separate_network_for_each_syllable
            try
                %% Plot the figure of errors for all networks over all trials.
                figure(9);
                replot_accuracies_concatanated('files', 'confusion_log_perf.txt');
            catch ME
            end
        end
        
        
        % Draw the hidden units' weights.  Let the user make these square or not
        % because lazy...
        if net.numLayers > 1 & false
            figure(5);
            for i = 1:size(net.IW{1}, 1)
                subplot(size(net.IW{1}, 1), 1, i)
                %imagesc([-time_window_steps:0]*fft_time_shift_seconds*1000, linspace(freq_range(1), freq_range(2), length(freq_range_ds))/1000, ...
                %    reshape(net.IW{1}(i,:), length(freq_range_ds), time_window_steps));
                imagesc([times(1:time_window_steps) - times(time_window_steps)]*1000, linspace(freq_range(1), freq_range(2), length(freq_range_ds))/1000, ...
                    reshape(net.IW{1}(i,:), length(freq_range_ds), time_window_steps));
                axis xy;
                ylabel('frequency');
                
                if i == 1
                    title('Hidden layers');
                end
                if i == size(net.IW{1}, 1)
                    xlabel('time (ms)');
                end
                %imagesc(reshape(net.IW{1}(i,:), time_window_steps, length(freq_range_ds)));
            end
        end
        drawnow;
        
        %save learn_detector_latest
        
        %% Save input file for the LabView detector
        % Extract data from net structure, because LabView's MathScript is too stupid to
        % permit the . operator.
        layer0 = net.IW{1};
        layer1 = net.LW{2,1};
        bias0 = net.b{1};
        bias1 = net.b{2};
        % The following store the transformation that took nnsetY to actual network training.  So
        % everything is forward.  Allow LabView to do:
        % out = (xmax-xmin)*(rawnetout-ymin)/(ymax-ymin) + xmin
        %     = (xrange/yrange) * (rawnetout - ymin) + xmin
        %     = (rawnetout - ymin) / gain + xmin
        mapstd_xmean = net.inputs{1}.processSettings{1}.xmean;
        mapstd_xstd = net.inputs{1}.processSettings{1}.xstd;
        if use_pattern_net
            mmmout_xmin = 0;
            mmmout_ymin = 0;
            mmmout_gain = 1;
        else
            mmmout_xmin = net.outputs{2}.processSettings{1}.xmin;
            mmmout_ymin = net.outputs{2}.processSettings{1}.ymin;
            mmmout_gain = net.outputs{2}.processSettings{1}.gain;
        end
        
        win_size = fft_size;
        fft_time_shift = fft_size - noverlap;
        scaling = 'linear';
        filename = sprintf('detector_%ss_frame%gms_%dhid_%dtrain.mat', ...
            sprintf('%g', toi), 1000*fft_time_shift_seconds_target, net.layers{1}.dimensions, ntrainsongsandnonsongs);
        fprintf('Saving as ''%s''...\n', filename);
        save(filename, ...
            'times_of_interest_s', 'toi', 'net', 'train_record', ...
            'samplerate', 'fft_size', 'win_size', 'fft_time_shift', 'fft_time_shift_seconds', 'fft_time_shift_seconds_target', ...
            'freq_range_ds', ...
            'time_window_steps', 'trigger_thresholds', 'freq_range', ...
            'layer0', 'layer1', 'bias0', 'bias1', ...
            'mmmout_xmin', 'mmmout_ymin', 'mmmout_gain', 'mapstd_xmean', 'mapstd_xstd', ...
            'shotgun_sigma', ...
            'ntrainsongsandnonsongs',  'scaling', '-v7');
        %% Save sample data: audio on channel0, canonical hits for first syllable on channel1
        if use_nn_realignment_test
            realignNetString = 'realignNet';
        else
            realignNetString = '';
        end
        
        % Let's standardise order for each set of test songs, so that we can compare multiple training
        % runs of the detector on the same songs:
        %rng(137);
        
        if create_song_test_file == 1 | (create_song_test_file == -1 & nruns == 1)
            if testfile_include_nonsinging
                % Re-permute all songs with a new random order
                newrand = randperm(size(mic_data,2));
                orig_songs_with_hits =  [ones(1, nmatchingsongs) zeros(1, nsongsandnonsongs - nmatchingsongs)]';
                new_songs_with_hits = orig_songs_with_hits(newrand);
                songs = reshape(mic_data(:, newrand), [], 1); % Include all singing and non-singing
                %songs = reshape(mic_data(:, 1:nsongsandnonsongs), [], 1); % Just singing
                songs_scale = max([max(songs) -min(songs)]);
                songs = songs / songs_scale;
                hits = zeros(size(mic_data));
                samples_of_interest = round(toi * samplerate);
                for i = 1:nsongsandnonsongs
                    if new_songs_with_hits(i)
                        % The baseline signal is recorded only for the first sample
                        % of interest:
                        hits(samples_of_interest(1) + sample_offsets_2(1, newrand(i)), i) = 1;
                    end
                end
                hits = reshape(hits, [], 1);
                songs = [songs hits];
                testfilename = sprintf('songs_%ss_%d%%%s.wav',...
                    sprintf('%g', toi(1)), round(100/(1+nonsinging_fraction)), ...
                    realignNetString);
            else
                % Just the real songs
                
                % Re-permute just 128 of the positive songs with a new random order -- for oscilloscope
                % 128-sample averages
                ntestsongsandnonsongs = nmatchingsongs;
                newrand = randperm(nmatchingsongs);
                newrand = newrand(1:ntestsongsandnonsongs);
                
                %songs = reshape(mic_data(:, 1:nmatchingsongs), [], 1); % Include all singing and non-singing
                songs = reshape(mic_data(:, newrand), [], 1); % Just singing
                songs_scale = max([max(songs) -min(songs)]);
                songs = songs / songs_scale;
                hits = zeros(nsamples_per_song, ntestsongsandnonsongs);
                samples_of_interest = round(toi * samplerate);
                for i = 1:ntestsongsandnonsongs
                    hits(samples_of_interest(1) + round(sample_offsets_test(1, i)), i) = 1;
                end
                hits = reshape(hits, [], 1);
                songs = [songs hits];
                
                testfilename = sprintf('songs_%d_%ss%s.wav', ...
                    ntestsongsandnonsongs, ...
                    sprintf('_%g', toi(1)), ...
                    realignNetString);
                
            end
            
            audiowrite(testfilename, songs, round(samplerate));
        end
        
        
        current_time = datetime('now');
        eta_date = start_time + (current_time - start_time) / ((run - first_run + 1) / (nruns - first_run + 1));
        if strcmp(datestr(eta_date, 'yyyymmdd'), datestr(current_time, 'yyyymmdd'))
            eta = datestr(eta_date, 15);
        else
            eta = datestr(eta_date, 'dddd yyyy-mm-dd HH:MM');
        end
        disp(sprintf('Expected finish time: %s.', eta));

    end
end

