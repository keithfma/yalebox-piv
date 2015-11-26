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
test_global_he = false;
test_local_he_brute_force = true;

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
    
    % plot original and equalized image
    figure
    
    subplot(2,1,1)
    imagesc(im);
    caxis([0,1]);
    axis off
    title('Original grayscale');
    
    subplot(2,1,2)
    imagesc(eql);
    caxis([0,1]);
    axis off
    title('Equalized');
    
    % plot empirical CDF and PDF for original and equalized images
    figure
    
    roi = im~=0;
    
    [cdf, val] = ecdf(im(roi));
    pdf = diff(cdf)./diff(val);
    
    subplot(2,2,1)
    plot(val, cdf, 'Marker', '.');
    title('Original CDF');
    
    subplot(2,2,3)
    plot(val(2:end), pdf, 'Marker', '.');
    title('Original PDF');
    
    [cdf, val] = ecdf(eql(roi));
    pdf = diff(cdf)./diff(val);
    
    subplot(2,2,2)
    plot(val, cdf, 'Marker', '.');
    title('Equalized CDF');
    
    subplot(2,2,4)
    plot(val(2:end), pdf, 'Marker', '.');
    title('Equalized PDF');
    
end

%% test 2: brute-force adaptive histogram equalization

if test_local_he_brute_force
    
    eql = local_he_brute_force(im, 0, brute_force_win);
    
end

% debug {
keyboard
% } debug
