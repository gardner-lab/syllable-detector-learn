% Place all params from the param file into an output struct. This lets us make sure they're all valid before importing them into
% the real workspace.
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
