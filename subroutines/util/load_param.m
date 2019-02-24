function param = load_param(param_file)
% function param = load_param(param_file)
% 
% Load parameters from JSON file
% %

param = load(param_file, '-mat');
