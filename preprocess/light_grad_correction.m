% function [correction] = light_grad_correction(img, mask, loess_frac);
% function [img_corrected] = light_grad_correction(img, mask, loess_frac);
%
% Compute lighting correction factor as a function of x for undesirable
% x-direction lighting gradients. The correction is applied as follows:
%
% img_corrected = img + repmat(correction, size(img,2), 1);
% img_corrected = min(img_max, max(img_min, img_corrected);
%
% where img_max and img_min define the valid range of image intensity. The
% correction is computed by (1) converting the image to grayscale, (2)
% extracting all the sand pixels using the mask, (3) fitting a smooth line
% to all extracted intensity data using loess, (4) computing the offset
% needed to flatten the smoothed trendline.
%
% Keith Ma, July 2015

% DEBUG: run as a script
img = imread('test_data/test_a_00_crop.png');
entropy_radius = 9;
open_radius = 36; 
bw_threshold = 0.75;
mask = sand_mask(img, entropy_radius, open_radius, bw_threshold);
loess_frac = 0.25;

%% check for sane inputs

assert(isa(img, 'uint8') && size(img,3) == 3, 'image is not 24-bit RGB');

assert(isa(mask, 'logical'), 'mask is not of type logical');
assert(size(mask,1) == size(img, 1) && size(mask,2) == size(img, 2), ...
    'mask and image dimensions do not match');

%% extract average intensity of sand pixels as (x, intensity) points

img = double(rgb2gray(img));
img = img-min(img(:));
img = img./max(img(:));
img(~mask) = NaN;
intens = nanmean(img,1);
 
%% fit smooth model

x = 1:size(img,2);
x_fit = downsample(x,10);
if x_fit(end) ~= x(end); x_fit(end+1) = x(end); end
 
intens_fit = loess(x, intens, x_fit, loess_frac, 1, 0);
 
% Note: higher-order interpolations differ by O(10^-6) from linear
intens_smooth = interp1(x_fit, intens_fit, x, 'linear');

%% compute correction function