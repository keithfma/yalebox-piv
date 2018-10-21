function [] = load_dependencies(varargin)
% [] = function load_dependencies(varargin)
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
    dep_name = varargin{ii};
    
    if strcmp(dep_name, 'jsonlab')
        add_dep(fullfile(base_path, 'depend', 'jsonlab-1.8'));
        
    elseif strcmp(dep_name, 'prep')
        add_dep(fullfile(base_path, 'sub_prep'));
        
        
    else
        error('Invalid dependency name: "%s"', dep_name);
        
    end
end

end
