function [] = post_series_from_file(post_file, piv_file, param_file)
% function [] = post_series_from_file(post_file, piv_file, param_file)
% 
% Post-processing for a PIV analysis output. Simple wrapper around
% post_series() that reads the input parameters from the standard JSON file
% rather than as arguments.
%
% Arguments:
%   post_file = String, filename of the MAT file to be created.
% 
%   piv_file = String, filename of the MAT file containing pre-processed
%       experiment images
%
%   param_file: JSON file containing all input parameters
% %

% load dependencies
update_path('jsonlab');

% load parameters struct
param = loadjson(param_file, 'SimplifyCell', 1);

% run post-processing with specified parameters
post_series(...
    piv_file, ...
    post_file, ...
    param.pro_bbox.value, ...
    param.retro_bbox.value, ...
    param.pad_method.value, ...
    param.notes.value);
