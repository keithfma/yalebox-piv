function [] = test_driver(prep)
% function [] = test_driver(prep)
%
% Run test for the experimental adaptive histogram equalization
% tools in this derectory. Most parameters are hard-coded in the first section.
%
% Arguments: 
%
% prep = Scalar, logical. Do image preparation steps (1) or load from existing
%          file (0).
% 
% %

%% initialize

% tests to run
test_global_he = 1;
test_local_he_brute_force = 0;

% environment
yalebox_piv_path = '/home/kfm/Documents/dissertation/yalebox-piv';
parpool_nworkers = 4;

% image prep
if nargin == 0; prep = true; end
input_file =  '/home/kfm/Documents/dissertation/yalebox-exp-fault/data/fault_ss_01/prep/2_crop/sidef/fault_ss_01_sidef_250.png'; %#ok!
output_file = 'test_image.mat'; % should be MAT
hue_lim = [0, 0.5]; 
val_lim = [0, 0.5];
entr_lim = [0.5, 1];
entr_win = 9;
morph_open_rad = 15;
morph_erode_rad = 10;

% equalization
brute_force_win = 31;

% start parallel pool if needed
if isempty(gcp('nocreate'))
    parpool(parpool_nworkers);
end

%% prepare image

addpath(yalebox_piv_path);
 
if prep == true
    
    % read image as rgb, hsv, and grayscale
    im_rgb = imread(image_file);
    im_hsv = rgb2hsv(im_rgb);
    im = im_hsv(:,:,3);
    
    % get and apply mask
    mask_manual = yalebox_prep_mask_manual(im_rgb);
    mask_auto = yalebox_prep_mask_auto(im_hsv, hue_lim, val_lim, entr_lim, ...
        entr_win, morph_open_rad, morph_erode_rad, true);
    im(~mask_auto | ~mask_manual) = 0;
    
    % save masked image for reuse, as .mat to preserve double values
    save(output_file, 'im');
    
else
    load(output_file, 'im');    
end

%% test 1: global histogram equalization

if test_global_he    
    eql = global_he(im, 0);
    display_test_results(im, 0, eql)
end

%% test 2: brute-force adaptive histogram equalization

if test_local_he_brute_force    
    eql = local_he_brute_force(im, 0, brute_force_win);
    display_test_results(im, 0, eql)    
end

% debug {
keyboard
% } debug 