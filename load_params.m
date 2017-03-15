% Place all params from the param file into an output struct. This lets us make sure they're all valid before importing them into
% the real workspace.
function [ p ] = load_params(params_file);
    
    if ~exist(strcat(params_file, '.m'), 'file')
        warning('Parameters file ''%s'' does not exist.', strcat(params_file, '.m'));
        p = {};
        return;
    end
    
    ignore = {'ans', 'v', 'params_file', 'ignore'};

    disp(sprintf('********** Configuration file %s: *************', ...
        strcat(pwd, filesep, params_file, '.m')));
    type(params_file);
    
    eval(params_file);
   
    v = whos;
    for i = 1:length(v)
        if isempty(intersect(v(i).name, ignore))
            eval(sprintf('p.%s = %s;', v(i).name, v(i).name));
        end
    end
end
