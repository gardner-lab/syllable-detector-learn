%% Plot the figure of errors for all networks over all trials...
function [] = replot_accuracies_concatanated(varargin)
    
    % The input file is created in show_confusion.m.  No effort is made to
    % ensure that it doesn't contain values for different configurations,
    % or even different-sized columns! So if you want to use it, best make
    % sure you start by deleting the previous confusion_log_perf.txt.  Keep
    % it in order to allow restart of partially completed jobs, since 10
    % syllables, 100 runs, 3000 training songs, convolutional networks,
    % etc., can take a long time to complete.
    
    real_bird_names = NaN;
    real_times = false;
    tex = false;
    
    for i = 1:2:nargin
        if ~exist(varargin{i}, 'var')
            warning('Argument ''%s'' is invalid.', varargin{i});
        else
            eval(sprintf('%s = varargin{i+1}', varargin{i}));
        end
    end
    
    % Custom code to pull out my datasets for publication:
    if ~exist('files', 'var') | isempty(files)
        files = {...
            '/Volumes/Data/song/lny64/confusion_log_perf_2ms.txt', ...
            '/Volumes/Data/song/lny64/confusion_log_perf.txt', ...
            '/Volumes/Data/song/lny42/confusion_log_perf_2ms.txt', ...
            '/Volumes/Data/song/lny42/confusion_log_perf.txt', ...
            '/Volumes/Data/song/lny46/confusion_log_perf_2ms.txt', ...
            '/Volumes/Data/song/lny46/confusion_log_perf.txt', ...
            '/Volumes/Data/song/lny4rb/confusion_log_perf_2ms.txt', ...
            '/Volumes/Data/song/lny4rb/confusion_log_perf.txt', ...
            '/Volumes/Data/song/lr12/all/confusion_log_perf_2ms.txt', ...
            '/Volumes/Data/song/lr12/all/confusion_log_perf.txt', ...
            '/Volumes/Data/song/lr13/all/confusion_log_perf_2ms.txt', ...
            '/Volumes/Data/song/lr13/all/confusion_log_perf.txt', ...
            '/Volumes/Data/song/lr28/all/confusion_log_perf_2ms.txt', ...
            '/Volumes/Data/song/lr28/all/confusion_log_perf.txt', ...
            '/Volumes/Data/song/lr77/all/confusion_log_perf_2ms.txt', ...
            '/Volumes/Data/song/lr77/all/confusion_log_perf.txt'
            };
        % If >= 0, this is the "/"-delimited field number for the name of the bird.
        real_bird_names = 4; % split() will count the leading "" before the first "/" as one, but I don't!
        real_times = true; % use actual detection times instead of t^*_n
    elseif ischar(files)
        foo = files;
        files = {foo};
    elseif iscell(files)
        % That's cool; do naught
    else
        error('The input parameter should be a confusion_log filename\nor a cell array thereof.');
    end
    
    if ~exist('tex', 'var')
        tex = false;
    end
    
    xtickl = {};
    xtime = {};
    xtabl = {};
    xfil = {};
    binj = [];
    confusion = [];
    
    n_so_far = 0;
    
    for f = 1:length(files)
        try
            confusions{f} = load(files{f});
        catch ME
            warning('Error loading ''%s''. Continuing...', files{f});
            continue;
        end
        
        % First, we have to do all this in order to count the unique syllables:
        [sylly bini binj_plus] = unique(confusions{f}(:,1));
        binj = [binj; binj_plus + n_so_far];
        
        for i = n_so_far+[1:length(sylly)]
            % Make the first column just a unique identifier, rather than the time target:
            id = 10+i;
            ids(i) = id;
            confusions{f}(find(confusions{f}(:,1) == sylly(i-n_so_far)), 1) = id;
            
            xtime{i} = sprintf('%g', 1000*sylly(i-n_so_far));
            % There's a table and a chart, bird names, bird autolabels ("A", "B", etc), times in ms, times in t^*_n...
            if i - n_so_far == 1
                if real_bird_names >= 0
                    foo = split(files{f}, {'\', '/'});
                    % split() will count the leading "" before the first "/" as one, but I don't! So add 1 here:
                    xtabl{i} = lower(strtok(foo(1+real_bird_names)));
                    if real_times
                        xtickl{i} = sprintf('%s_{%s}', xtabl{i}, xtime{i});
                    else
                        xtickl{i} = sprintf('%s t^*_%d', xtabl{i}, i-n_so_far);
                    end
                else
                    if real_times
                        xtickl{i} = sprintf('%c_{%s}', 'A'+f-1, xtime{i});
                    else
                        xtickl{i} = sprintf('%c t^*_%d', 'A'+f-1, i-n_so_far);
                    end
                    xtabl{i} = sprintf('%c', 'A'+f-1);
                end
                xfil{i} = files{f};
            else
                if real_times
                    xtickl{i} = xtime{i};
                else
                    xtickl{i} = sprintf('t^*_%d', i-n_so_far);
                end
                xtabl{i} = ' ';
                xfil{i} = '';
            end
            sylly_counts(i) = length(find(confusions{f}(:,1)==id));
            sylly_means(i,:) = mean(confusions{f}(find(confusions{f}(:,1)==id),2:3), 1);
        end
        
        n_so_far = n_so_far + length(sylly);
        confusion = [confusion ; confusions{f}];
    end
    
    performance = sprintf('Bird\tTime\tTrue Pos %%\tFalse Pos %%\tRuns\tFile');
    for i = 1:length(xtickl)
        if tex
            if i == length(xtickl)
                lineterm = '';
            else
                lineterm = '\\';
            end
            str = sprintf('%s & \t%s & \t%s & \t%s %s %% %s', xtabl{i}, xtime{i}, ...
                sigfig(sylly_means(i,1)*100, 4), ...
                sigfig(sylly_means(i,2)*100, 2, 'pad'), ...
                lineterm, ...
                xfil{i});
            performance(i+1, 1:length(str)) = str;
        else
            str = sprintf('%s\t%s\t%s\t\t%s\t\t%d\t%s', xtabl{i}, xtime{i}, ...
                sigfig(sylly_means(i,1)*100, 4), ...
                sigfig(sylly_means(i,2)*100, 2, 'pad'), ...
                sylly_counts(i), ...
                xfil{i});
            performance(i+1, 1:length(str)) = str;
        end
    end
    performance
    
    colours = distinguishable_colors(length(xtickl));
    offsets = (rand(size(confusion(:,1))) - 0.5) * 2 * 0.33;
    if size(confusion, 2) >= 4 & false
        sizes = (mapminmax(-confusion(:,4)')'+1.1)*8;
    else
        sizes = 3;
    end
    try
        subplot(1,2,1);
        scatter(confusion(:,1)+offsets, confusion(:,2)*100, sizes, colours(binj,:), 'filled');
        xlabel('Test syllable');
        ylabel('True Positives %');
        title('Correct detections');
        %if min(confusion(:,1)) ~= max(confusion(:,1))
        %    set(gca, 'xlim', [min(confusion(:,1))-0.4 max(confusion(:,1))+0.4]);
        %end
        %set(gca, 'ylim', [97 100]);
        set(gca, 'xtick', ids, 'xticklabel', xtickl);
        if exist('xticklabel_rotate', 'file')
            xticklabel_rotate([], 60);
        end
        
        subplot(1,2,2);
        scatter(confusion(:,1)+offsets, confusion(:,3)*100, sizes, colours(binj,:), 'filled');
        xlabel('Test syllable');
        ylabel('False Positives %');
        title('Incorrect detections');
        %if length(ids) > 1
        %    set(gca, 'xlim', [min(ids)-0.5 max(ids)+0.5]);
        %end
        set(gca, 'xtick', ids, 'xticklabel', xtickl);
        if exist('xticklabel_rotate', 'file')
            xticklabel_rotate([], 60);
        end
    end
    %set(gca, 'ylim', [0 0.07]);
    %sylly_counts
end
