% Script. Template script for exploring and selecting image pre-processing
% parameters. The intended wusage is to make a copy of the script for a given
% experiment, run it cell by cell, modifying the default parameters to suit the
% experiment particulars.

%% Init

woco_image_file = 'test/prep_get_parameters_woco.jpg';
step_image_file = 'test/prep_get_parameters_step.jpg';
param_out_file = 'test/prep_get_parameters.mat';

%% Define coordinate system control points

[ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp] = ...
    prep_world_coord_control_pts(woco_image_file, true);

%% Rectify and crop 

% parameters
crop_xw = [-0.600, 0.600];
crop_yw = [ 0.002, 0.200];

% actions
[woco, xw, yw] = ...
    prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, ...
        imread(woco_image_file), true);

[rgb, ~, ~] = ...
    prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, ...
        imread(step_image_file), true);

%% Define a manual mask

mask_manual = prep_mask_manual(rgb);

%% Apply automatic masking model

hue_lim = [0.01, 0.3];
val_lim = [0.05, 0.5];
entropy_lim = [0.4, 1];
entropy_len = 11;
morph_open_rad = 10;
morph_erode_rad = 5;

mask_auto = prep_mask_auto(rgb, hue_lim, val_lim, entropy_lim, entropy_len, ...
                    morph_open_rad, morph_erode_rad, true, true);

%% adaptive histogram equalization

% parameters
eql_len = 21;

% actions
eql = prep_intensity(rgb, mask_manual & mask_auto, eql_len, true, true);

%% Get image file directory and file list

output_file = 'junk.image.nc';
image_path = '/home/kfm/Documents/dissertation/yalebox-exp-erosion/data/k24/image/clean';
image_name_glob = 'K24*.jpg';

tmp = dir(fullfile(image_path, image_name_glob));
image_names = {tmp(:).name};

%% Save parameters for batch processing

save(param_out_file, 'ctrl_xw', 'ctrl_yw', 'ctrl_xp', 'ctrl_yp', 'crop_xw', ...
    'crop_yw', 'entropy_len', 'num_cluster', 'cluster_center', 'eql_len', ...
    'xw', 'yw', 'mask_manual', 'output_file', 'image_path', 'image_names'); 
