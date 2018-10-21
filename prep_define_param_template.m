% Define parameters for image pre-processing
%
% Usage: run this script to populate the workspace with the defined parameters,
% then run the prep_param.m script to test the results and save them when you
% are satified. 
%
% This file is a TEMPLATE -- make a copy and edit the values for your
% particular experiment details.
% %


% image files parameters --------------------------------------------------
%
% woco_image_file: string, path to world coordinate image
% step_image_file: string, path to step image to use for testing
% %
woco_image_file = '/home/keith/prj/yalebox-faultless-wedge/fault_ss_01_sidef_clean/fault_ss_01_woco_sidef.jpg';
step_image_file = '/home/keith/prj/yalebox-faultless-wedge/fault_ss_01_sidef_clean/fault_ss_01_sidef_320.jpg';


% world coordinate control points parameters ------------------------------
%
% ctrl_pts_action: string, select action to perform for control points,
%   choose from 'try', 'retry', 'load'
% ctrl_pts_mat_file: string, path to MAT file to load/save world coordinate
%   control points (depending on selected action)
% %
ctrl_pts_mat_file = '/home/keith/prj/yalebox-faultless-wedge/fault_ss_01_sidef_prep/fault_ss_01_sidef_prep_ctrl_pts.mat';
ctrl_pts_action = 'try'; 


% rectify and crop --------------------------------------------------------
%
% crop_xw: 2-element array, minimum and maximum x coordinates for crop in
%   world coordinate units (m)
% crop_yw: 2-element array, minimum and maximum y coordinates for crop in
%   world coordinate units (m)
% %
crop_xw = [-0.5, 0.65];
crop_yw = [ 0.002, 0.150];


% manual image masking ----------------------------------------------------
%
% mask_manual_action: string, select action to perform for manual mask,
%   choose from 'create', 'load'
% mask_manual_mat_file: string, path to MAT file to load/save manual mask
% %
mask_manual_action = 'create';
mask_manual_file = '/home/keith/prj/yalebox-faultless-wedge/fault_ss_01_sidef_prep/fault_ss_01_sidef_prep_mask_manual.mat';


% automatic image masking -------------------------------------------------
%
% hue_lim: See help for prep_mask_auto() 
% value_lim: " "
% entropy_lim: " "
% entropy_len: " "
% morph_open_rad: " "
% morph_erode_rad: " "
% %
hue_lim = [0.02, 0.2];
value_lim = [0.15, 1.0];
entropy_lim = [0.6, 1];
entropy_len = 11;
morph_open_rad = 10;
morph_erode_rad = 2;
