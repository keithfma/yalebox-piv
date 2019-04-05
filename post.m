function [] = post(param, piv_mat, post_mat)
% function [] = post(param, piv_mat, post_mat)
% 
% Post-processing for a PIV analysis output. Simple wrapper around
% post_series() that reads the input parameters from the standard JSON file
% rather than as arguments.
%
% Arguments:
% 
%   param: Either (a) string, path to JSON file containing all input
%       parameters, or (b) struct containing same
%
%   piv_mat = String, filename of the MAT file containing pre-processed
%       experiment images
%
%   post_mat = String, filename of the MAT file to be created.
% %

update_path('util');

switch class(param)
    case 'struct'
        param = param;  % do nothing
    case 'char'
        param = load_param(param);
    otherwise
        error('Unexpected type for input argument "param"');
end

post_series(...
    piv_mat, ...
    post_mat, ...
    param.pro_bbox.value, ...
    param.retro_bbox.value, ...
    param.pad_method.value, ...
    param.notes.value);
