function [] = save_param(param, param_file)
% function [] = save_param(param, param_file)
% 
% Write parameters to JSON file
% 
% Arguments:
%   param: parameter struct to be saved
%   param_file: string, path to save parameters as JSON
% %

save(param_file, '-mat', '-v7.3', '-struct', 'param');