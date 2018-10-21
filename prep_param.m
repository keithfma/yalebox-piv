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


%% Define coordinate system control points ---------------------------------
% Create or edit world coordinate control points interactively. You will be
% prompted to save the results by updating the param JSON file upon
% completion.
% %

param = loadjson(PREP_PARAM_FILE);
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
    fprintf('"woco" parameters NOT written to file');
end


%%  Rectify and crop 
% Apply rectification using control points, and crop specified by
% parameters, then display the results
% % 

param = loadjson(PREP_PARAM_FILE);

[woco, xw, yw] = ...
    prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, ...
        imread(woco_image_file), true);

[rgb, ~, ~] = ...
    prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, ...
        imread(step_image_file), true);

%% Manual image mask 
% Either create or load manual image mask
% % 

if strcmp(mask_manual_action, 'create')
    % define mask interactively
    mask_manual = prep_mask_manual(rgb);
    % decide whether to save, task care to avoid accidental overwriting
    save_mask_manual = true;
    if exist(mask_manual_mat_file, 'file') ~= 2
        button = questdlg('Manual mask file exists! Overwrite it?', 'WARNING', 'Yes', 'No', 'No');
        if strcmp(button, 'No'); save_mask_manual = false; end
    end
    % save mask to MAT file
    if save_mask_manual
        save(mask_manual_mat_file, 'mask_manual');
    end

elseif strcmp(mask_manual_action, 'load')
    % load mask from previous results
    load(mask_manual_mat_file)

else
    error('Invalid choice for parameter "mask_manual_action"');

end

% % Automatic image mask ---------------------------------------------------
% % apply automatic image masking routine given specified parameters
% % % 
% 
% mask_auto = prep_mask_auto(rgb, hue_lim, value_lim, entropy_lim, entropy_len, ...
%                     morph_open_rad, morph_erode_rad, true, true);


%% TODO; run prep for several images in the series