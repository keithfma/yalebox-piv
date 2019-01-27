% Define and/or test pre-processing parameters
%
% Expects path to the prep parameters definition file defined in a variable
% PREP_PARAM_FILE. To get started, call prep_default_param to create a 
% template JSON file with default parameter values, then edit the variables
% therein to suit your experiment.
%
% This script is designed to be run cell-by-cell. Each cell runs one step
% of the image preparation process. Run a cell, inspect the results, and
% edit the parameter file until you are satisfied with the results. Note
% that some cells depend on the results of previous cells.
% %

update_path('prep', 'jsonlab', 'util');

%% define "images" section parameters -- interactive

param = load_param(PREP_PARAM_FILE);

% path to images directory
path = uigetdir(...
    param.images.path.value, 'Select directory containing image files');
path = strip(path, 'right', filesep);

% world coordinate image
[woco_file, woco_path] = uigetfile('*.*', 'Select world coordinate image', ...
    fullfile(param.images.path.value, param.images.woco_file.value));
woco_path = strip(woco_path, 'right', filesep);
assert(strcmp(woco_path, path), 'Expect image to be in image directory');

% test image to use for parameter definition in this script
[test_file, test_path] = uigetfile('*.*', 'Select test image', ...
    fullfile(param.images.path.value, param.images.test_file.value));
test_path = strip(test_path, 'right', filesep);
assert(strcmp(test_path, path), 'Expect image to be in selected image directory');

% all experiment images
[exp_files, exp_path] = uigetfile(...
    '*.*', 'Select all experiment images', path, 'MultiSelect', 'on');
exp_path = strip(exp_path, 'right', filesep);
assert(strcmp(exp_path, path), 'Expect images to be in selected image directory');

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write parameters to %s?', PREP_PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.images.path.value = path;
    param.images.woco_file.value = woco_file;
    param.test.file.value = test_file;
    param.images.exp_files.value = exp_files;
    save_param(param, PREP_PARAM_FILE);
    fprintf('Saved parameter file: %s\n', PREP_PARAM_FILE);
else
    fprintf('Parameter file NOT saved\n');
end

%% define "woco" section parameters -- interactive

param = load_param(PREP_PARAM_FILE);

% sanity check
lengths = [length(param.woco.xw.value), ...
    length(param.woco.yw.value), ...
    length(param.woco.xp.value), ...
    length(param.woco.yp.value)];
if min(lengths) ~= max(lengths)
    error('Length of control point vectors do not all match');
end

% define control points, either from scratch or starting from existing
if max(lengths) > 0
    [xw, yw, xp, yp, ~] = prep_world_coord_control_points(...
        fullfile(param.images.path.value, param.images.woco_file.value), ...
        param.woco.xw.value, ...
        param.woco.yw.value, ...
        param.woco.xp.value, ...
        param.woco.yp.value, ...
        true(size(param.woco.xw.value)));
else
    [xw, yw, xp, yp, ~] = prep_world_coord_control_points(...
        fullfile(param.images.path.value, param.images.woco_file.value));
end

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write parameters to %s?', PREP_PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.woco.xw.value = xw;
    param.woco.yw.value = yw;
    param.woco.xp.value = xp;
    param.woco.yp.value = yp;
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
        param.woco.xp.value, ...
        param.woco.yp.value, ...
        param.woco.xw.value, ...
        param.woco.yw.value, ...
        param.crop.xlim.value, ...
        param.crop.ylim.value, ...
        imread(fullfile(param.images.path.value, param.images.woco_file.value)), ...
        [], ...
        true);

fprintf('Displaying rectifed & cropped test image\n');
[rgb, ~, ~] = ...
    prep_rectify_and_crop(...
        param.woco.xp.value, ...
        param.woco.yp.value, ...
        param.woco.xw.value, ...
        param.woco.yw.value, ...
        param.crop.xlim.value, ...
        param.crop.ylim.value, ...
        imread(fullfile(param.images.path.value, param.test.file.value)), ...,
        [], ...
        true);

%% define "mask_manual" section parameters -- interactive

param = load_param(PREP_PARAM_FILE);

% interactive mask creation
[mask_manual, poly] = prep_mask_manual(...
    rgb, ...
    param.mask_manual.poly.value, ...
    true);

% create masked image for next step
mrgb = rgb;
mrgb(repmat(~mask_manual, 1, 1, 3)) = 0;

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write parameters to %s?', PREP_PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.mask_manual.poly.value = poly;
    save_param(param, PREP_PARAM_FILE);
    fprintf('Saved parameter file: %s\n', PREP_PARAM_FILE);
else
    fprintf('Parameter file NOT saved\n');
end

%% apply "mask_auto" section parameters -- edit param file to update

param = load_param(PREP_PARAM_FILE);

mask_auto = prep_mask_auto(...
    mrgb, ...
    param.mask_auto.hue_lim.value,...
    param.mask_auto.value_lim.value, ...
    param.mask_auto.entropy_lim.value, ...
    param.mask_auto.entropy_len.value, ...
    param.mask_auto.morph_open_rad.value, ...
    param.mask_auto.morph_erode_rad.value, ...
    true, ...
    true);

%% apply histogram equalization

param = load_param(PREP_PARAM_FILE);

prep_intensity(...
    rgb, ...
    mask_manual & mask_auto, ...
    param.intensity.eql_len.value, ...
    true, ...
    true);

%% run small-scale test case

param = load_param(PREP_PARAM_FILE);

% select evenly-spaced files from exp file list
idx = round(linspace(...
    1, ...
    length(param.images.exp_files.value), ...
    param.test.num_images.value));
param.images.exp_files.value = param.images.exp_files.value(idx);

% save temporary param file, register for cleanup on function completion
temp_param_file = get_temp_file('json');
clean_param = onCleanup(@() delete(temp_param_file));
savejson('', param, 'Filename', temp_param_file, 'SingletArray', 0);

% call prep series
result_file = get_temp_file('nc');
fprintf('Test run results: %s\n', result_file);
prep_series_from_file(result_file, temp_param_file);

% plot results for each image in sequence
test_img = ncread(result_file, 'img');
test_img_rgb = ncread(result_file, 'img_rgb');
test_mask_auto = ncread(result_file, 'mask_auto');
test_mask_manual = ncread(result_file, 'mask_manual');

hf = figure;
for ii = 1:param.test.num_images.value
    hf.Name = sprintf('Test Analysis Results: %s', ...
        param.images.exp_files.value{ii});
    subplot(3,1,1)
    imshow(test_img_rgb(:,:,:,ii));
    subplot(3,1,2)
    imshow(test_mask_manual & test_mask_auto(:,:,ii))
    subplot(3,1,3)
    imshow(test_img(:,:,ii));
    pause
end
 