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
%   color_tol = Scalar.
%   smooth_rad = Scalar.    
%
% Keith Ma, July 2015

% DEBUG: Run as a script with predefined inputs
image_file = 'test_data/test_a_00.jpg';
color_ref = [];
color_tol = [];
smooth_rad = [];
