function [] = update_path(varargin)
% [] = function update_path(varargin)
% 
% Load specified dependencies to the matlab path
%
% Arguments:
%
% Returns: Nothing
% % 

[base_path, ~, ~] = fileparts(mfilename('fullpath'));

function [] = add_dep(dep_path)
    addpath(dep_path);
    fprintf('Added %s to MATLAB path\n', dep_path);
end

for ii = 1:length(varargin)
    
    switch varargin{ii}
        
        % third-party libraries
        case 'jsonlab'
            add_dep(fullfile(base_path, 'dependencies', 'jsonlab-1.8'));
    
        % subroutines
        case 'prep'
            add_dep(fullfile(base_path, 'subroutines', 'prep'));
        case 'piv'
            add_dep(fullfile(base_path, 'subroutines', 'piv'));
        case 'util'
            add_dep(fullfile(base_path, 'subroutines', 'util'));
        
        % parameter definition
        case 'check_prep'
            add_dep(fullfile(base_path, 'subroutines', 'check_prep'));
        
        % bad name, fail
        otherwise
            error('Invalid dependency name: "%s"', dep_name);
        
    end
end

end