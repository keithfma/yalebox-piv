% Utility script to help in setting / assessing image prep parameters
%
% This script is designed to be run cell-by-cell. Run a cell, inspect the
% results, and edit the parameter file until you are satisfied with the
% results. Note that some cells may depend on the results of previous
% cells.
%
% Expected variables in workspace:
% 
%   PARAM_FILE: path to parameters definition file. Use
%       piv_default_param() to create a template, and populate the
%       variables therein to suit your experiment.
% 
%   TEST_INDEX: int, index of initial image to use for single step test
% %

update_path('prep', 'util');

fprintf('Running prep parameter check with:\n');
fprintf('\tPARAM_FILE = %s\n', PARAM_FILE);
fprintf('\tTEST_INDEX = %i\n', TEST_INDEX);

%% define "images" section parameters -- interactive

param = load_param(PARAM_FILE);

% path to images directory
path = uigetdir(...
    param.image_dir.value, 'Select directory containing image files');
path = strip(path, 'right', filesep);

% world coordinate image
[woco_file, woco_path] = uigetfile('*.*', 'Select world coordinate image', ...
    fullfile(param.image_dir.value, param.woco_file.value));
woco_path = strip(woco_path, 'right', filesep);
assert(strcmp(woco_path, path), 'Expect image to be in image directory');

% mask training image
[train_file, train_path] = uigetfile('*.*', 'Select mask training image', ...
    fullfile(param.image_dir.value, param.mask_train_file.value));
train_path = strip(train_path, 'right', filesep);
assert(strcmp(train_path, path), 'Expect image to be in image directory');

% all experiment images
[exp_files, exp_path] = uigetfile(...
    '*.*', 'Select all experiment images', path, 'MultiSelect', 'on');
exp_path = strip(exp_path, 'right', filesep);
assert(strcmp(exp_path, path), 'Expect images to be in selected image directory');

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write parameters to %s?', PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.image_dir.value = path;
    param.woco_file.value = woco_file;
    param.mask_train_file.value = train_file;
    param.exp_files.value = exp_files;
    save_param(param, PARAM_FILE);
    fprintf('Saved parameter file: %s\n', PARAM_FILE);
else
    fprintf('Parameter file NOT saved\n');
end

%% define "woco" section parameters -- interactive

param = load_param(PARAM_FILE);

% sanity check
lengths = [length(param.woco_xw.value), ...
    length(param.woco_yw.value), ...
    length(param.woco_xp.value), ...
    length(param.woco_yp.value)];
if min(lengths) ~= max(lengths)
    error('Length of control point vectors do not all match');
end

% define control points, either from scratch or starting from existing
if max(lengths) > 0
    [xw, yw, xp, yp, ~] = prep_world_coord_control_points(...
        fullfile(param.image_dir.value, param.woco_file.value), ...
        param.woco_xw.value, ...
        param.woco_yw.value, ...
        param.woco_xp.value, ...
        param.woco_yp.value, ...
        true(size(param.woco_xw.value)));
else
    [xw, yw, xp, yp, ~] = prep_world_coord_control_points(...
        fullfile(param.image_dir.value, param.images.woco_file.value));
end

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write parameters to %s?', PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.woco_xw.value = xw;
    param.woco_yw.value = yw;
    param.woco_xp.value = xp;
    param.woco_yp.value = yp;
    save_param(param, PARAM_FILE);
    fprintf('Saved parameter file: %s\n', PARAM_FILE);
else
    fprintf('Parameter file NOT saved\n');
end

%% apply "crop" parameters -- edit param file to update

param = load_param(PARAM_FILE);

fprintf('Displaying rectifed & cropped woco image\n');
[~, ~, ~] = ...
    prep_rectify_and_crop(...
        param.woco_xp.value, ...
        param.woco_yp.value, ...
        param.woco_xw.value, ...
        param.woco_yw.value, ...
        param.crop_xlim.value, ...
        param.crop_ylim.value, ...
        imread(fullfile(param.image_dir.value, param.woco_file.value)), ...
        true);

fprintf('Displaying rectifed & cropped test image\n');
[rgb, ~, ~] = ...
    prep_rectify_and_crop(...
        param.woco_xp.value, ...
        param.woco_yp.value, ...
        param.woco_xw.value, ...
        param.woco_yw.value, ...
        param.crop_xlim.value, ...
        param.crop_ylim.value, ...
        imread(fullfile(param.image_dir.value, param.exp_files.value{TEST_INDEX})), ...,
        true);

%% define labels for mask training -- interactive

param = load_param(PARAM_FILE);

[train_rgb, ~, ~] = prep_rectify_and_crop(...
        param.woco_xp.value, ...
        param.woco_yp.value, ...
        param.woco_xw.value, ...
        param.woco_yw.value, ...
        param.crop_xlim.value, ...
        param.crop_ylim.value, ...
        imread(fullfile(param.image_dir.value, param.mask_train_file.value)));

fprintf('Label training data for all classes\n');
[~, mask_poly] = prep_mask_labels(train_rgb, param.mask_train_labels.value, true);

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write parameters to %s?', PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.mask_train_labels.value = mask_poly;
    save_param(param, PARAM_FILE);
    fprintf('Saved parameter file: %s\n', PARAM_FILE);
else
    fprintf('Parameter file NOT saved\n');
end

%% extract features for mask training -- edit param file to update

% TODO: specify parameters in parameter file
scale = 25;
sigma = 0.5;
min_size = 50;

train_segments = prep_mask_segments(train_rgb, scale, sigma, min_size, true);

train_features = prep_mask_features(train_rgb, train_segments);


%% train and apply sand/other masking model  -- edit param file to update

param = load_param(PARAM_FILE);

mask_model = prep_mask_train(...
    train_rgb, ...
    param.mask_train_labels.value, ...
    param.mask_model_type.value);

mask = prep_mask_apply(rgb, mask_model, true);

    
%% define "mask_manual" section parameters -- interactive

param = load_param(PARAM_FILE);

% interactive mask creation
[mask_manual, poly] = prep_mask_manual(...
    rgb, ...
    param.mask_poly.value, ...
    true);

% create masked image for next step
mrgb = rgb;
mrgb(repmat(~mask_manual, 1, 1, 3)) = 0;

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write parameters to %s?', PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.mask_poly.value = poly;
    save_param(param, PARAM_FILE);
    fprintf('Saved parameter file: %s\n', PARAM_FILE);
else
    fprintf('Parameter file NOT saved\n');
end

%% apply "mask_auto" section parameters -- edit param file to update

param = load_param(PARAM_FILE);

mask_auto = prep_mask_auto(...
    mrgb, ...
    param.mask_hue_lim.value,...
    param.mask_value_lim.value, ...
    param.mask_entropy_lim.value, ...
    param.mask_entropy_len.value, ...
    param.mask_morph_open_rad.value, ...
    param.mask_morph_erode_rad.value, ...
    true);

%% apply histogram equalization

param = load_param(PARAM_FILE);

prep_intensity(...
    rgb, ...
    mask_manual & mask_auto, ...
    param.intensity_eql_len.value, ...
    true);

 