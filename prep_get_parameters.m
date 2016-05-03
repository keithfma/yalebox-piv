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
[rgb, xw, yw] = ...
    prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, ...
        imread(step_image_file), true);

%% Define a manual mask

mask_manual = prep_mask_manual(rgb);

%% Train and apply automatic masking model

% parameters
entropy_len = 11;
num_cluster = 3;

% actions
cluster_center = prep_mask_auto_train(rgb, entropy_len, num_cluster, true, true);

%% Masked adaptive histogram equalization

% parameters
eql_len = 21;

% actions
mask_auto = prep_mask_auto(rgb, entropy_len, cluster_center, true, true);
hsv = rgb2hsv(rgb);
value = hsv(:,:,3);
eql = prep_intensity(value, mask_manual & mask_auto, eql_len, true, true);

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