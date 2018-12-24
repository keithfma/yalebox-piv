function [] = prep_series_from_file(result_file, param_file)
% function [] = prep_series_from_file(result_file, param_file)
% 
% Create PIV input file for a given image series. Reads in the images,
% rectifies and crops, masks, corrects illumination, and saves the results
% and metadata in a netCDF file.
%
% Simple wrapper around prep_series() that reads the (many) input
% parameters from the standard JSON file rather than as arguments.
%
% Arguments:
%   result_file = String, filename of the netCDF input file to be created. 
%   param_file: JSON file containing all input parameters
% %

% load dependencies
update_path('jsonlab', 'prep');

% load parameters struct
param = loadjson(param_file, 'SimplifyCell', 1);

% generate missing input parameters
% note: not changed in prep_series in order to maintain backward compat
img = imread(...
    fullfile(param.images.path.value, param.images.exp_files.value{1}));
% missing world coordinate vectors
[~, xw, yw] = prep_rectify_and_crop(...
    param.woco.xp.value, ...
    param.woco.yp.value, ...
    param.woco.xw.value, ...
    param.woco.yw.value, ...
    param.crop.xlim.value, ...
    param.crop.ylim.value, ...
    img, ...
    param.crop.npts.value, ...
    false, ...
    false);
% missing manual mask
[mask_manual, ~] = prep_mask_manual(...
    zeros(length(yw), length(xw), 3), ...
    param.mask_manual.poly.value, ...
    false);


% call standard prep_series function
prep_series(...
    result_file, ...
    param.images.path.value, ...
    param.images.exp_files.value, ...
    param.woco.xw.value, ...
    param.woco.yw.value, ...
    param.woco.xp.value, ...
    param.woco.yp.value, ...
    param.crop.xlim.value, ...
    param.crop.ylim.value, ...
    param.crop.npts.value, ...
    param.mask_auto.hue_lim.value, ...
    param.mask_auto.value_lim.value, ...
    param.mask_auto.entropy_lim.value, ...
    param.mask_auto.entropy_len.value, ...
    param.mask_auto.morph_open_rad.value, ...
    param.mask_auto.morph_erode_rad.value, ...
    param.intensity.eql_len.value, ...
    xw, ...
    yw, ...
    mask_manual);
