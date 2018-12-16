% Test pre-processing parameters and save to MAT file if satisfactory
%
% Expects path to the prep parameters definition file defined in a variable
% PREP_PARAM_FILE. Use the template prep_param_template.json as a starting
% point, and populate the variables therein to suit your experiment.
%
% This script is designed to be run cell-by-cell. Each cell runs one step
% of the image preparation process. Run a cell, inspect the results, and
% edit the parameter file until you are satisfied with the results. Note
% that some cells depend on the results of previous cells.
% %

%% Load dependencies
% Add required dependencies to the MATLAB path, you MUST run this cell once
% for the subsequent cells to work
% %

load_dependencies('prep', 'jsonlab');


%% Define paths to experiment images

param = loadjson(PREP_PARAM_FILE, 'SimplifyCell', 1);

% world coordinate image
woco_file = uigetfile(...
    '*.*', 'Select World Coordinate Image', param.images.woco_file.value);

% test image to use for parameter definition in this script
test_file = uigetfile(...
    '*.*', 'Select Test Image', param.images.woco_file.value);

% all experiment images
exp_files = uigetfile(...
    '*.*', 'Select All Experiment Images', 'MultiSelect', 'on');

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write "images" parameters to %s?', PREP_PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.images.woco_file.value = woco_file;
    param.images.test_file.value = test_file;
    param.images.exp_files.value = exp_files;
    savejson('', param, 'Filename', PREP_PARAM_FILE, 'SingletArray', 0);
    fprintf('"images" parameters written to file: %s\n', PREP_PARAM_FILE);
else
    fprintf('"images" parameters NOT written to fileopt.SimplifyCell\n');
end



%% Define coordinate system control points ---------------------------------
% Create or edit world coordinate control points interactively. You will be
% prompted to save the results by updating the param JSON file upon
% completion.
% %

param = loadjson(PREP_PARAM_FILE, 'SimplifyCell', 1);
woco = param.woco;
woco_file = param.images.woco_file.value;

% check for existing control points, check sanity
if ~isempty(woco.xw) || ~isempty(woco.yw) || ~isempty(woco.xp) || ~isempty(woco.yp)
    npts = length(woco.xw.value);
    if length(woco.yw.value) ~= npts || length(woco.xp.value) ~= npts || length(woco.yp.value) ~= npts
        error('Length of control point vectods do not match');
    end
    retry = true;
else
    retry = false;
end
    
% define control points, either from scratch or starting from existing
if retry
    [xw, yw, xp, yp, ~] = prep_world_coord_control_points(...
        woco_file, woco.xw.value, woco.yw.value, woco.xp.value, ...
        woco.yp.value, true(size(woco.xw.value)));
else
    [xw, yw, xp, yp, ~] = prep_world_coord_control_points(woco_file);
end
    
% save results, take care to avoid accidental overwriting
prompt = sprintf('Write "woco" parameters to %s?', PREP_PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.woco.xw.value = xw;
    param.woco.yw.value = yw;
    param.woco.xp.value = xp;
    param.woco.yp.value = yp;
    savejson('', param, 'Filename', PREP_PARAM_FILE, 'SingletArray', 0);
    fprintf('"woco" parameters written to file: %s\n', PREP_PARAM_FILE);
else
    fprintf('"woco" parameters NOT written to fileopt.SimplifyCell\n');
end


%%  Rectify and crop 
% Apply rectification using control points, and crop specified by
% parameters, then display the results
% % 

param = loadjson(PREP_PARAM_FILE, 'SimplifyCell', 1);
images = param.images;
woco = param.woco;
crop = param.crop;

fprintf('Displaying rectifed & cropped woco image using current parameters\n');
[~, xw, yw] = ...
    prep_rectify_and_crop(...
        woco.xp.value, woco.yp.value, woco.xw.value, woco.yw.value, ...
        crop.xlim.value, crop.ylim.value, imread(images.woco_file.value), ...
        [], true);

fprintf('Displaying rectifed & cropped test image using current parameters\n');
[rgb, ~, ~] = ...
    prep_rectify_and_crop(...
        woco.xp.value, woco.yp.value, woco.xw.value, woco.yw.value, ...
        crop.xlim.value, crop.ylim.value, imread(images.test_file.value), ...,
        [], true);
    

%% Manual image mask 
% Interactively edit manual image mask
% %

% load parameters
param = loadjson(PREP_PARAM_FILE, 'SimplifyCell', 1);

% interactive mask creation
[mask_manual, poly] = prep_mask_manual(...
    rgb, param.mask_manual.poly.value, true);

% create masked image for next step
mrgb = rgb;
mrgb(repmat(~mask_manual, 1, 1, 3)) = 0;

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write "mask_manual" parameters to %s?', PREP_PARAM_FILE);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.mask_manual.poly.value = poly;
    savejson('', param, 'Filename', PREP_PARAM_FILE, 'SingletArray', 0);
    fprintf('"mask_manual" parameters written to file: %s\n', PREP_PARAM_FILE);
else
    fprintf('"mask_manual" parameters NOT written to file\n');
end


%% Automatic image mask 
% apply automatic image masking routine given specified parameters
% % 

% load parameters
param = loadjson(PREP_PARAM_FILE, 'SimplifyCell', 1);
mask = param.mask_auto;

% apply auto masking and display results
mask_auto = prep_mask_auto(...
    rgb, mask.hue_lim.value, mask.value_lim.value, ...
    mask.entropy_lim.value, mask.entropy_len.value, ...
    mask.morph_open_rad.value, mask.morph_erode_rad.value, ...
    true, true);


%% TODO; run prep for several images in the series
