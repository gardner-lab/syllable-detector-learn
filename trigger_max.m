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
function [ trigger_now val pos] = trigger_max(responses, threshold, varargin);

% De-bounce the signal

% This just returns the best hit per song...
responses(find(responses<threshold)) = NaN;

[ val pos ] = max(responses, [], 2);

trigger_now = zeros(size(responses));

for i = 1:length(pos)
    if isnan(val(i))
        pos(i) = NaN;
    else
        trigger_now(i, pos(i)) = 1;
    end
end
