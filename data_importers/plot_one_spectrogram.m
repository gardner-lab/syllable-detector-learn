% Plot a spectrogram from microphone data.
%   mic_data:  audio data
%   fs:        sampling frequency
%   threshold: OPTIONAL detection threshold; if included, also add 'scores', 'starts' and 'ends'
%   scores:    Dynamic Timewarping match scores
%   starts:    points at which ostensible matches start (seconds)
%   ends:      points at which ostensible matches end (seconds)

function plot_one_spectrogram(mic_data, fs, threshold, scores, starts, ends);

    
    
    fft_size = 256;                                  % FFT size
    fft_time_shift_seconds_target = 0.004;
    
    nsamples_per_song = length(mic_data);
    song_seconds = nsamples_per_song / fs;
    sample_times = 0:1/fs:song_seconds;
    sample_times = sample_times(1:nsamples_per_song);
    
    noverlap = fft_size - (floor(fs * fft_time_shift_seconds_target));
    % SPECGRAM(A,NFFT=512,Fs=[],WINDOW=[],noverlap=500)
    %speck = specgram(mic_data(:,1), 512, [], [], 500) + eps;
    
    window = hamming(fft_size);
    
    [speck freqs times] = spectrogram(mic_data, window, noverlap, [], fs);
    [nfreqs, ntimes] = size(speck);
    speck = speck + eps;
    
    % This will be approximately the same as fft_time_shift_seconds_target, but not quite: the fft_time_shift
    % is given by noverlap, and will actually be fft_size/target_samplerate
    speck = abs(speck);
    
    %% Draw the pretty full-res spectrogram and the targets
    figure(412);
    imagesc([times(1) times(end)]*1000, [freqs(1) freqs(end)]/1000, log(speck));
    set(gca, 'YLim', [0 10]);
    axis xy;
    xlabel('Time (ms)');
    ylabel('Frequency (kHz)');
    colormap bone;
    
    if exist('threshold', 'var')
        for i = 1:length(starts)
            if scores(i) > threshold
                a=line([starts(i) ends(i)]*1000, [0.5 0.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 10);
            else
                a=line([starts(i) ends(i)]*1000, [0.5 0.5], 'Color', [0 0.8 0.4], 'LineWidth', 20);
            end
        end
    end
    
    drawnow;
    %disp('Press any key...');
    %pause;
end
