%function [mask] = create_mask(image_file, color_ref, color_tol, smooth_rad)
% function [mask] = create_mask(image_file, color_ref, color_tol, smooth_rad)
%
% Create a logical mask for a color image (image_file) that is TRUE where
% there is sand and FALSE elsewhere. This can be used to remove (set to 0)
% the background in a image prior to PIV analysis or other applications.
%
% Sand is identified by a reference color (color_ref, the average color or
% sand) in (RGB? CMYK? Other?) colorspace, a tolerance (color_tol, a
% normalized distance to the reference color in ND colorspace), and a
% smoothing radius (smooth_rad, radius of moving average filter in pixels).
%
% A pixel is sand if:
%
% ||smoothed_pixel(r,g,b) - color_ref(r,g,b)|| < color_tol
%
% Arguments:
%
%   image_file = String. Filename of image to be analyzed.
%   color_ref = 3 element vector.
%   color_tol = Scalar. Color tolerance as a normalized distance in ()
%               colorspace, valid range is [0,1]
%   smooth_rad = Scalar. Radius of smoothing filter in pixels.    
%
% Keith Ma, July 2015

% DEBUG: Run as a script with predefined inputs
image_file = 'test_data/test_a_00.jpg';
color_ref = [0, 0, 0];
color_tol = 0.5;
smooth_rad = 10;

%% check for sane inputs

assert(exist(image_file, 'file') == 2, 'image_file is not a valid filename');

tmp = size(color_ref);
assert(numel(tmp) == 2 && min(tmp) == 1 && max(tmp) ==3, ...
    'color_ref is not a 3-element vector');

assert(numel(color_tol) == 1, 'color_tol is not a scalar');
assert(color_tol >= 0 && color_tol <= 1, 'color_tol is not in the range [0,1]');

assert(numel(smooth_rad) == 1, 'smooth_rad is not a scalar');

%% read in data