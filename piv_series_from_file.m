function [] = piv_series_from_file(result_file, image_file, param_file)
% function [] = piv_series_from_file(result_file, param_file)
% 
% Run PIV analysis for a given image series. Simple wrapper around
% piv_series() that reads the (many) input parameters from the standard
% JSON file rather than as arguments.
%
% Arguments:
%   result_file = String, filename of the MAT file to be created.
% 
%   image_file = String, filename of the MAT file containing pre-processed
%       experiment images
%
%   param_file: JSON file containing all input parameters
% %

% load dependencies
update_path('jsonlab');

% load parameters struct
param = loadjson(param_file, 'SimplifyCell', 1);

% run PIV with specified parameters
piv_series(...
    result_file, ...
    image_file, ...
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
