% Utility script to help in setting / assessing image prep parameters
%
% This script is designed to be run cell-by-cell. Run a cell, inspect the
% results, and edit the parameter file until you are satisfied with the
% results. Note that some cells may depend on the results of previous
% cells.
%
% Expected variables in workspace:
% 
%   PREP_PARAM_FILE: path to parameters definition file. Use
%       piv_default_param() to create a template, and populate the
%       variables therein to suit your experiment.
% 
%   TEST_INDEX: int, index of initial image to use for single step test
% %

update_path('prep', 'util');

fprintf('Running prep parameter check with:\n');
fprintf('\tPREP_PARAM_FILE = %s\n', PREP_PARAM_FILE);
fprintf('\tTEST_INDEX = %i\n', TEST_INDEX);

%% define "images" section parameters -- interactive

param = load_param(PREP_PARAM_FILE);

% path to images directory
path = uigetdir(...
    param.image_dir.value, 'Select directory containing image files');
path = strip(path, 'right', filesep);

% world coordinate image
[woco_file, woco_path] = uigetfile('*.*', 'Select world coordinate image', ...
    fullfile(param.image_dir.value, param.woco_file.value));
woco_path = strip(woco_path, 'right', filesep);
assert(strcmp(woco_path, path), 'Expect image to be in image directory');

% all experiment images
[exp_files, exp_path] = uigetfile(...
    '*.*', 'Select all experiment images', path, 'MultiSelect', 'on');
exp_path = strip(exp_path, 'right', filesep);
assert(strcmp(exp_path, path), 'Expect images to be in selected image directory');

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write parameters to %s?', PREP_PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.image_dir.value = path;
    param.woco_file.value = woco_file;
    param.exp_files.value = exp_files;
    save_param(param, PREP_PARAM_FILE);
    fprintf('Saved parameter file: %s\n', PREP_PARAM_FILE);
else
    fprintf('Parameter file NOT saved\n');
end

%% define "woco" section parameters -- interactive

param = load_param(PREP_PARAM_FILE);

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
prompt = sprintf('Write parameters to %s?', PREP_PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.woco_xw.value = xw;
    param.woco_yw.value = yw;
    param.woco_xp.value = xp;
    param.woco_yp.value = yp;
    save_param(param, PREP_PARAM_FILE);
    fprintf('Saved parameter file: %s\n', PREP_PARAM_FILE);
else
    fprintf('Parameter file NOT saved\n');
end

%% apply "crop" parameters -- edit param file to update

param = load_param(PREP_PARAM_FILE);

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
[rgb, x, y] = ...
    prep_rectify_and_crop(...
        param.woco_xp.value, ...
        param.woco_yp.value, ...
        param.woco_xw.value, ...
        param.woco_yw.value, ...
        param.crop_xlim.value, ...
        param.crop_ylim.value, ...
        imread(fullfile(param.image_dir.value, param.exp_files.value{TEST_INDEX})), ...,
        true);

    
%% define "mask_manual" section parameters -- interactive

param = load_param(PREP_PARAM_FILE);

% interactive mask creation
[mask_manual, poly] = prep_mask_manual(...
    rgb, ...
    param.mask_poly.value, ...
    true);

% create masked image for next step
mrgb = rgb;
mrgb(repmat(~mask_manual, 1, 1, 3)) = 0;

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write parameters to %s?', PREP_PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.mask_poly.value = poly;
    save_param(param, PREP_PARAM_FILE);
    fprintf('Saved parameter file: %s\n', PREP_PARAM_FILE);
else
    fprintf('Parameter file NOT saved\n');
end

%% apply "mask_auto" section parameters -- edit param file to update

param = load_param(PREP_PARAM_FILE);

mask_auto = prep_mask_auto(...
    mrgb, ...
    param.mask_hue_lim.value,...
    param.mask_value_lim.value, ...
    param.mask_entropy_lim.value, ...
    param.mask_entropy_len.value, ...
    param.view.value, ...
    true);

%% apply "equalize" section parameters -- edit param file to update

tic;

param = load_param(PREP_PARAM_FILE);

img_eql = prep_equalize(...
    rgb, ...
    mask_auto & mask_manual, ...
    param.equalize_len.value, ...
    true);

toc

%% apply pad and fill -- edit param file to update

param = load_param(PREP_PARAM_FILE);

[x_pad, y_pad, img_pad, mask_pad] = prep_pad(...
    x, ...
    y, ...
    img_eql, ...
    mask_auto & mask_manual, ...
    param.pad_num_rows.value, ...
    param.pad_num_cols.value);

img_fill = prep_fill(...
    img_pad, ...
    mask_pad, ...
    param.fill_method.value, ...
    true);
