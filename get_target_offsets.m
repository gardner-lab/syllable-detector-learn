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
function [ target_offsets ] = get_target_offsets(MIC_DATA, times_of_interest, fs, timestep_length_ds, canonical_songs);


% How big a window is 20ms?
window_halfsize = floor(0.03 * fs);

nsongs = size(MIC_DATA, 2);

samples_of_interest = floor(times_of_interest * timestep_length_ds * fs);

target_offsets = zeros(length(times_of_interest), nsongs);

for i = 1:length(samples_of_interest)
        template = MIC_DATA(samples_of_interest(i) - window_halfsize : samples_of_interest(i) + window_halfsize, ...
                            canonical_songs(i));
        for j = 1:nsongs
                cor = xcorr(template, ...
                            MIC_DATA(samples_of_interest(i) - window_halfsize : samples_of_interest(i) + window_halfsize, j));
                [val pos] = max(cor);
                target_offsets(i,j) = pos - 2 * window_halfsize - 1;
        end
end

% Convert from sample steps to timestep_length_ds units
target_offsets = round( ( target_offsets / fs ) / timestep_length_ds );
