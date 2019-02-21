% Utility script to help in setting / assessing PIV post-processing parameters
%
% This script is designed to be run cell-by-cell. Run a cell, inspect the
% results, and edit the parameter file until you are satisfied with the
% results. Note that some cells may depend on the results of previous
% cells.
%
% Expected variables in workspace:
% 
%   PARAM_FILE: path to parameters definition file. Use
%       post_default_param() to create a template, and populate the
%       variables therein to suit your experiment.
% 
%   PIV_FILE: path to MAT-file output from PIV analysis
%
%   TEST_INDEX: int, index of initial image to use for single step test
% %

update_path('post', 'jsonlab', 'util');

fprintf('Running post parameter check with:\n');
fprintf('\tPARAM_FILE = %s\n', PARAM_FILE);
fprintf('\tPIV_FILE = %s\n', PIV_FILE);
fprintf('\tTEST_INDEX = %i\n', TEST_INDEX);


%% pro- and retro- sampling boxes

param = load_param(PARAM_FILE);
piv = matfile(PIV_FILE, 'Writable', false);

xx = piv.x_grd;
yy = piv.y_grd;
mask = double(piv.roi_grd(:, :, TEST_INDEX));
mask(~mask) = NaN;
uu = piv.u_grd(:, :, TEST_INDEX).*mask;
vv = piv.v_grd(:, :, TEST_INDEX).*mask;

post_displ_rect(xx, yy, uu, vv, param.pro_bbox.value, true);
hf = gcf;
hf.Name = 'Pro';

post_displ_rect(xx, yy, uu, vv, param.retro_bbox.value, true);
hf = gcf;
hf.Name = 'Retro';
