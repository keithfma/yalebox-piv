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
%   result_file = String, filename of the MAT input file to be created. 
%   param_file: JSON file containing all input parameters
% %

update_path('jsonlab', 'prep');

pp = loadjson(param_file, 'SimplifyCell', 1);

prep_series(...
    result_file, ...
    pp.image_dir.value, ...
    pp.exp_files.value, ...
    pp.view.value, ...
    pp.woco_xw.value, ...
    pp.woco_yw.value, ...
    pp.woco_xp.value, ...
    pp.woco_yp.value, ...
    pp.crop_xlim.value, ...
    pp.crop_ylim.value, ...
    pp.crop_npts.value, ...
    pp.mask_poly.value, ...
    pp.mask_hue_lim.value, ...
    pp.mask_value_lim.value, ...
    pp.mask_entropy_lim.value, ...
    pp.mask_entropy_len.value, ...
    pp.mask_morph_open_rad.value, ...
    pp.mask_morph_erode_rad.value, ...
    pp.intensity_eql_len.value, ...
    pp.notes.value);
