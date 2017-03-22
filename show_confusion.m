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
function [ ] = show_confusion(...
        testout, ...
        nwindows_per_song, ...
        FALSE_POSITIVE_COST, ...
        times_of_interest, ...
        tstep_of_interest, ...
        MATCH_PLUSMINUS, ...
        timestep, ...
        time_window_steps, ...
        songs_with_hits, ...
        trigger_thresholds);


% A positive that happens within ACTIVE_TIME of the event does not count as a
% false positive.  This is in seconds, and allows for some jitter.
ACTIVE_TIMESTEPS_BEFORE = floor(MATCH_PLUSMINUS / timestep);
ACTIVE_TIMESTEPS_AFTER = floor(MATCH_PLUSMINUS / timestep);

% The timesteps of interest are with reference to the start of the song.
% Responses have been trimmed to start at the start of recognition given
% the time window.  So we need to align those:
tstep_of_interest_shifted = tstep_of_interest - time_window_steps + 1;

for i = 1:length(tstep_of_interest)
    responses = squeeze(testout(i, :, :))';
    positive_interval = tstep_of_interest_shifted(i)-ACTIVE_TIMESTEPS_BEFORE:...
        tstep_of_interest_shifted(i)+ACTIVE_TIMESTEPS_AFTER;
    positive_interval = positive_interval(find(positive_interval > 0 & positive_interval <= nwindows_per_song));
    
    f = @(threshold)trigger_threshold_cost(threshold, ...
        responses, ...
        positive_interval, ...
        FALSE_POSITIVE_COST, ...
        songs_with_hits);
    
    % Find optimal threshold on the interval [0.001 1]
    % optimal_thresholds = fminbnd(f, 0.001, 1);
    % fminbnd is useless at jumping out of local minima (and this will have large flat regions), but brute-forcing the search is quick.
    [ ~, trueposrate, falseposrate ] = f(trigger_thresholds(i));

    tpfp = [times_of_interest trueposrate falseposrate 0];

    if true
        % ASCII art
        fprintf('At %d ms:        True positive    negative\n', times_of_interest(i) * 1000);
        fprintf('     output pos      %.5f%%     %s%%\n', trueposrate*100, sigfig(falseposrate*100));
        fprintf('            neg       %s%%       %.5f%%\n', sigfig((1-trueposrate)*100), (1-falseposrate)*100);
        if exist('confusion_log_perf.txt', 'file')
            save('confusion_log_perf.txt', 'tpfp', '-append', '-ascii');
        else
            save('confusion_log_perf.txt', 'tpfp', '-ascii');
        end
    else
        % LaTeX table
        fprintf('\\vspace{8pt}\\par\\noindent\n\\begin{tabular}{r|cc}\n  {\\bf At %d ms} & \\multicolumn{2}{c}{True} \\\\ \n  & pos & neg \\\\ \n  \\hline  Detected pos & %.5f\\%% & %.5f\\%%\\\\ \n  neg & %.5f\\%% & %.5f\\%%\\\\ \n\\end{tabular}\n', ...
            times_of_interest(i) * 1000, trueposrate*100, falseposrate*100, (1-trueposrate)*100, (1-falseposrate)*100);
    end
end

