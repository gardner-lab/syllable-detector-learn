%% Plot the figure of errors for all networks over all trials...
figure(9);
% This file is created in show_confusion.m.  No effort is made to ensure that it doesn't
% contain values for different configurations, or even different-sized columns! So if
% you want to use it, best make sure you start by deleting the previous
% confusion_log_perf.txt.  Keep it in order to allow restart of partially completed
% jobs, since 10 syllables, 100 runs, 3000 training songs, convolutional networks, etc.,
% can take a long time to complete.

files = {...%'~/v/syllable-detector-learn/data/confusion_log_perf_4hid.txt', ...
    '/Volumes/Data/song/lny64/confusion_log_perf.txt', ...
    '/Volumes/Data/song/LNY46 0.28s or 0.195s/confusion_log_perf.txt', ...
    '/Volumes/Data/song/LNY42 0.38s/confusion_log_perf.txt', ...
    ...%'/Volumes/Data/song/LNY4RB 0.25s/confusion_log_perf_400nonsong.txt', ...
    '/Volumes/Data/song/LNY4RB 0.25s/confusion_log_perf.txt', ...
    '/Volumes/Data/song/lr28/all/confusion_log_perf.txt', ...
    '/Volumes/Data/song/lr13/all/confusion_log_perf.txt', ...
    '/Volumes/Data/song/lr12/all/confusion_log_perf.txt' ...
    };


xtickl = {};
sylly_means = [];
sylly_counts = [];
confusion = [];
binj = [];

n_so_far = 0;

for f = 1:length(files)
    confusions{f} = load(files{f});
    
    % First, we have to do all this in order to count the unique syllables:
    [sylly bini binj_plus] = unique(confusions{f}(:,1));
    binj = [binj; binj_plus + n_so_far];
    
    for i = n_so_far+[1:length(sylly)]
        % Make the first column just a unique identifier, rather than the time target:
        id = 10+i;
        ids(i) = id;
        confusions{f}(find(confusions{f}(:,1) == sylly(i-n_so_far)), 1) = id;
        
        xtickl{i} = sprintf('%c t^*_%d', 'A'+f-1, i-n_so_far);
        sylly_counts(i) = length(find(confusions{f}(:,1)==id));
        sylly_means(i,:) = mean(confusions{f}(find(confusions{f}(:,1)==id),2:3), 1);
    end
    
    n_so_far = n_so_far + length(sylly);
    confusion = [confusion ; confusions{f}];
end

sylly_means_percent = sylly_means * 100
colours = distinguishable_colors(length(xtickl));
offsets = (rand(size(confusion(:,1))) - 0.5) * 2 * 0.33;
if size(confusion, 2) >= 4 & false
    sizes = (mapminmax(-confusion(:,4)')'+1.1)*8;
else
    sizes = 3;
end
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

subplot(1,2,2);
scatter(confusion(:,1)+offsets, confusion(:,3)*100, sizes, colours(binj,:), 'filled');
xlabel('Test syllable');
ylabel('False Positives %');
title('Incorrect detections');
%if length(ids) > 1
%    set(gca, 'xlim', [min(ids)-0.5 max(ids)+0.5]);
%end
set(gca, 'xtick', ids, 'xticklabel', xtickl);
%set(gca, 'ylim', [0 0.07]);
sylly_counts
