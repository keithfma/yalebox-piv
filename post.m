function [] = post(param_json, piv_mat, post_mat)
% function [] = post(param_json, piv_mat, post_mat)
% 
% Post-processing for a PIV analysis output. Simple wrapper around
% post_series() that reads the input parameters from the standard JSON file
% rather than as arguments.
%
% Arguments:
% 
%   param_json: JSON file containing all input parameters
% 
%   piv_mat = String, filename of the MAT file containing pre-processed
%       experiment images
%
%   post_mat = String, filename of the MAT file to be created.
% %

update_path('jsonlab');
param = loadjson(param_json, 'SimplifyCell', 1);
post_series(...
    piv_mat, ...
    post_mat, ...
    param.pro_bbox.value, ...
    param.retro_bbox.value, ...
    param.pad_method.value, ...
    param.notes.value);
