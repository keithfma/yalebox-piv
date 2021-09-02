function new_experiment(folder, varargin)
% function new_experiment(folder, varargin) 
%
% Populate a new directory with all the files you need to start processing a new experiment.
%
% Arguments:
%   folder: the destination folder
%
% Optional name-value Parameters:
%   mkdir: set true to create the folder if it does not already exist (default false)
%   overwrite: set true to overwrite existing files (default false)
% % 

p = inputParser();
p.FunctionName = 'new_experiment';
addRequired(p, 'folder', @(x) validateattributes(x, {'char'}, {'vector'}));
addParameter(p, 'mkdir', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
addParameter(p, 'overwrite', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
parse(p, folder, varargin{:});

fprintf('New experiment in folder: %s (mkdir: %d, overwrite: %d)\n', ...
    folder, p.Results.mkdir, p.Results.overwrite);

% root of the yalebox-piv software package
YALEBOX_ROOT_PATH = fileparts(mfilename('fullpath'));

% template file paths, to be copied as live scripts (.mlx), assume each of the template
%   file paths end with _template.m
TEMPLATE_FILES = {...
    fullfile(YALEBOX_ROOT_PATH, 'subroutines', 'prep', 'prep_param_template.m') ...
};

% create directory, if requested and required
if p.Results.mkdir && ~isfolder(folder)
    mkdir(folder);
end


% copy template files and covert to live-scripts, overwritting if requested and required
for idx = 1:length(TEMPLATE_FILES)
    
    template_file = TEMPLATE_FILES{idx};
    
    [~, template_name, ~] = fileparts(template_file); 
    dest_file = fullfile(folder, [template_name,  '.mlx']); 
    
    if ~p.Results.overwrite && isfile(dest_file)
        error('Destination file exists and overwriting is disabled: %s', dest_file);
    end
    
    fprintf('\t%s\n', dest_file);
    matlab.internal.liveeditor.openAndSave(template_file, dest_file);
    
end
