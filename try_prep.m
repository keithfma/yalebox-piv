% Define and/or test pre-processing parameters
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

load_dependencies('prep', 'prep_param', 'jsonlab', 'util');

%% define "images" section parameters -- interactive

define_images(PREP_PARAM_FILE);

%% define "woco" section parameters -- interactive

define_woco(PREP_PARAM_FILE);

%% apply "crop" parameters -- edit param file to update

rgb = apply_crop(PREP_PARAM_FILE);

%% define "mask_manual" section parameters -- interactive

[mrgb, mask_manual] = define_mask_manual(PREP_PARAM_FILE, rgb);

%% apply "mask_auto" section parameters -- edit param file to update

mask_auto = apply_mask_auto(PREP_PARAM_FILE, mrgb);

%% apply histogram equalization

apply_intensity(PREP_PARAM_FILE, mrgb, mask_manual & mask_auto);

%% run small-scale test case

test_run(PREP_PARAM_FILE);
 