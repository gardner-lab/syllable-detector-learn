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
function [ tstep_of_interest ] = suggest_moments_of_interest(ntsteps, ...
                                             spectrogram_avg_img, ...
                                             time_window_steps, ...
                                             timestep, ...
                                             layer0sz, ...
                                             nwindows_per_song, ...
                                             ntimes_ds, ...
                                             freq_range_ds);
                                             

if 0  % Look for unique moments
        dists = squareform(pdist(spectrogram_avg_img(:,time_window_ds:end)'));
else   % Look for unique segments from the training set
        pdistset = zeros(layer0sz, nwindows_per_song);
        for tstep = time_window_steps : ntimes_ds
                pdistset(:, tstep - time_window_steps + 1) ...
                        = reshape(spectrogram_avg_img(freq_range_ds, ...
                        tstep - time_window_steps + 1  :  tstep), ...
                        [], 1);
        end

        dists = squareform(pdist(pdistset'));
end

figure(6);
imagesc(dists);

min_sep_s = 0.05;
min_sep_samples = ceil(timestep / min_sep_s)

%% We want to find the row such that the minimum distance to that point is
%% maximised--not counting points close to the row.  So set the diagonal high.
for i = 1:size(dists, 1)
        for j = i:size(dists, 1)
                if abs(i-j) < min_sep_samples
                        dists(i, j) = NaN;
                        dists(j, i) = NaN;
                end                        
        end
end

% Cut off the ends
dists(:,1:min_sep_samples) = NaN;
dists(:,end-min_sep_samples+1:end) = NaN;
dists(1:min_sep_samples,:) = NaN;
dists(end-min_sep_samples+1:end) = NaN;


tstep_of_interest = [];

%% Pick the maximally-distal spot on the distance image, add it to the
%% timesteps we might care about, delete it from the image, repeat.
for i = 1:ntsteps
        [val pos] = max(min(dists));
        tstep_of_interest(i) = pos;
        
        dists(pos-min_sep_samples:pos+min_sep_samples,:) = NaN * zeros(2*min_sep_samples+1, nwindows_per_song);
        dists(:, pos-min_sep_samples:pos+min_sep_samples) = NaN * zeros(nwindows_per_song, 2*min_sep_samples+1);
        hold on;
        scatter(pos, pos, 100, 'r*');
        hold off;
end
tstep_of_interest = tstep_of_interest + time_window_steps - 1
tstep_of_interest = sort(tstep_of_interest);

drawnow;
