function [] = piv(param, prep_mat, piv_mat)
% function [] = piv(param, prep_mat, piv_mat)
% 
% Run PIV analysis for a given image series. Simple wrapper around
% piv_series() that reads the (many) input parameters from the standard
% MAT file rather than as arguments.
%
% Arguments:
% 
%   param: Either (a) string, path to MAT file containing all input
%       parameters, or (b) struct containing same
% 
%   prep_mat = String, filename of the MAT file containing pre-processed
%       experiment images
% 
%   piv_mat = String, filename of the MAT file to be created.
% %

update_path('util', 'piv');

switch class(param)
    case 'struct'
        param = param;  %#ok! do nothing
    case 'char'
        param = load_param(param);
    otherwise
        error('Unexpected type for input argument "param"');
end
        
piv_series(...
    piv_mat, ...
    prep_mat, ...
    param.step_range.value, ...
    param.gap.value, ...
    param.samp_len.value, ...
    param.samp_spc.value, ...
    param.intr_len.value, ...
    param.num_pass.value, ...
    param.valid_radius.value, ...
    param.valid_max.value, ...
    param.valid_eps.value, ...
    param.min_frac_data.value, ...
    param.min_frac_overlap.value, ...
    param.notes.value);