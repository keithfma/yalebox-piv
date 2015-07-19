%function [mask] = create_mask(image_file, color_ref, color_tol, smooth_rad)
% function [mask] = create_mask(image_file, color_ref, color_tol, smooth_rad)
%
% Create a logical mask for a color image (image_file) that is TRUE where
% there is sand and FALSE elsewhere. This can be used to remove (set to 0)
% the background in a image prior to PIV analysis or other applications.
%
% Sand is identified by a reference color (color_ref, the average color or
% sand) in 24-bit RGB (e.g. 24-bit Truecolor), a tolerance (color_tol, a
% maximum distance to the reference color in normalized RGB colorspace),
% and a smoothing radius (smooth_rad, radius of moving average filter in
% pixels).
%
% A pixel is sand if:
%
% ||smoothed_pixel(r,g,b) - color_ref(r,g,b)|| / 255 < color_tol
%
% Arguments:
%
%   image_file = String. Filename of image to be analyzed.
%   color_ref = 3 element vector. 24-bit RGB representation of the average
%               color of sand, valid range is 0-255.
%   color_tol = Scalar. Color tolerance as distance in normalized RGB
%               colorspace (where 0-255 is rescaled to 0-1 for each axis), valid range is [XXX]
%   smooth_rad = Scalar. Radius of smoothing filter in pixels, can be fractional.    
%
% Keith Ma, July 2015

% DEBUG: Run as a script with predefined inputs
image_file = 'test_data/test_a_00_crop.png';
color_ref = [28.9640   28.6109   24.7218];
color_tol = 0.5;
smooth_rad = 10;

%% read in data, check for sane inputs

assert(exist(image_file, 'file') == 2, 'image_file is not a valid filename');
img = imread(image_file);
assert(isa(img, 'uint8') && size(img,3) == 3, 'image is not 24-bit RGB');

tmp = size(color_ref);
assert(numel(tmp) == 2 && min(tmp) == 1 && max(tmp) == 3, ...
    'color_ref is not a 3-element vector');
assert(max(color_ref) <= 255, min(color_ref) >= 0, ...
    'color_ref is not in the range [0. 255]');

assert(numel(color_tol) == 1, 'color_tol is not a scalar');
assert(color_tol >= 0 && color_tol <= 1, 'color_tol is not in the range [0,1]');

assert(numel(smooth_rad) == 1, 'smooth_rad is not a scalar');


%% smooth image

%filt = fspecial('disk', smooth_rad);
%img = imfilter(img, filt);

N = 20;
img = double(img)/255;
for i = 1:3
    img(:,:,i) = medfilt2(img(:,:,i));
end

% %% normalize colors
% 
% img = double(img)/255;
% color_ref = color_ref/255;
% 
% %% compute distances to reference color
% 
% img_dist = sqrt( (img(:,:,1)-color_ref(1)).^2 + (img(:,:,2)-color_ref(2)).^2 + (img(:,:,3)-color_ref(3)).^2 );







