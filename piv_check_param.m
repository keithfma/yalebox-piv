% Test PIV parameters
%
% Expected variables:
%   PIV_PARAM_FILE: path to parameters definition file. Use
%       templates/piv.json as a starting point, and populate the variables
%       therein to suit your experiment.
%   IMAGES_FILE: path to pre-processed image file netCDF, as produced by
%       prep_series.m
%
% This script is designed to be run cell-by-cell. Each cell runs one step
% of the image preparation process. Run a cell, inspect the results, and
% edit the parameter file until you are satisfied with the results. Note
% that some cells depend on the results of previous cells.
% %

update_path('jsonlab', 'xcorr', 'spline', 'deriv', 'piv', 'post', 'util');

%% PIV analysis for single image pair -- edit parameter file to change

param = load_param(PIV_PARAM_FILE);

% read image data
[xw, yw, ini, ini_mask, fin, fin_mask] = read_image_pair_from_nc(...
    IMAGES_FILE, ...
    param.test.ini.value, ...
    param.piv.gap.value);

% compute PIV
piv_result = piv(...
    ini, ...
    fin, ...
    ini_mask, ...
    fin_mask, ...
    xw, ...
    yw, ...
    param.piv.samp_len.value, ...
    param.piv.samp_spc.value, ...
    param.piv.intr_len.value, ...
    param.piv.num_pass.value, ...
    param.piv.valid_radius.value, ...
    param.piv.valid_max.value, ...
    param.piv.valid_eps.value, ...
    param.piv.min_frac_data.value, ...
    param.piv.min_frac_overlap.value, ...
    true);

% compute strain (helps to judge PIV result quality)
strain_result = post_strain(...
    piv_result.x_grd(1,:), ...
    piv_result.y_grd(:,1), ... 
    piv_result.u_grd, ...
    piv_result.v_grd, ...
    piv_result.roi_grd, ...
    'nearest'); 

% display results for user validation
display_image_pair(xw, yw, ini, ini_mask, fin, fin_mask);
display_piv(piv_result);
display_strain(strain_result);

%% run small-scale test case

param = load_param(PIV_PARAM_FILE);

% TODO: start here!

% % select evenly-spaced files from exp file list
% idx = round(linspace(...
%     1, ...
%     length(param.images.exp_files.value), ...
%     param.test.num_images.value));
% param.images.exp_files.value = param.images.exp_files.value(idx);
% 
% % save temporary param file, register for cleanup on function completion
% temp_param_file = get_temp_file('json');
% clean_param = onCleanup(@() delete(temp_param_file));
% savejson('', param, 'Filename', temp_param_file, 'SingletArray', 0);
% 
% % call prep series
% result_file = get_temp_file('nc');
% fprintf('Test run results: %s\n', result_file);
% prep_series_from_file(result_file, temp_param_file);
% 
% % plot results for each image in sequence
% test_img = ncread(result_file, 'img');
% test_img_rgb = ncread(result_file, 'img_rgb');
% test_mask_auto = ncread(result_file, 'mask_auto');
% test_mask_manual = ncread(result_file, 'mask_manual');
% 
% hf = figure;
% for ii = 1:param.test.num_images.value
%     hf.Name = sprintf('Test Analysis Results: %s', ...
%         param.images.exp_files.value{ii});
%     subplot(3,1,1)
%     imshow(test_img_rgb(:,:,:,ii));
%     subplot(3,1,2)
%     imshow(test_mask_manual & test_mask_auto(:,:,ii))
%     subplot(3,1,3)
%     imshow(test_img(:,:,ii));
%     pause
end