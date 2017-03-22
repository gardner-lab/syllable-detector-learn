% Place all params from the param file into an output struct. This lets us make sure they're all valid before importing them into
% the real workspace.

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
function [ p ] = load_params(params_file);
    
    if ~exist(strcat(params_file, '.m'), 'file')
        warning('Parameters file ''%s'' does not exist.', strcat(params_file, '.m'));
        p = {};
        return;
    end
    
    % The following are temporary variables created by this function. Don't save them.
    ignore = {'ans', 'v', 'params_file', 'ignore'};
    
    eval(params_file);
   
    v = whos;
    for i = 1:length(v)
        if isempty(intersect(v(i).name, ignore))
            eval(sprintf('p.%s = %s;', v(i).name, v(i).name));
        end
    end
end
