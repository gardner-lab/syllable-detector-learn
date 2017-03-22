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
function [ trigger_now ] = trigger(responses, threshold, schmidt_trigger_down, timestep);

% De-bounce the signal


% Look at all trigger events:
trigger_now = responses > threshold;
trigger_schmidt = zeros(size(responses));

% How about a Schmidt Trigger?
if exist('schmidt_trigger_down')
        schmidt_trigger_down_steps = schmidt_trigger_down / timestep;
        for i = 1:size(responses, 1)
                pos = find(trigger_now(i,:) > 0);
                while ~isempty(pos)
                        trigger_schmidt(i, min(pos)) = 1;
                        pos(find(pos <= min(pos) + schmidt_trigger_down_steps)) = [];
                end
        end
        trigger_now = trigger_schmidt;
end

