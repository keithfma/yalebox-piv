%function [color_ref, subset] = get_color_ref(image_file, subset)
% function [color_ref, subset] = get_color_ref(image_file, subset)
%
%
% Return the mean color (24-bit RGB) of a subset of the the image in
% image_file. The subset region can be selected interactively (only
% image_file is provided) or can be speficied in the function call (both
% image_file and subset are provided).
%
% Arguments:
%
%   image_file = String. Filename of the image to be analyzed.
%   subset = 4-element vector. Limits of the subset to be averaged, in the
%           form [min_row, max_row, min_col, max_col]. Will be set
%           interactively if not provided as an input argument.
%   color_ref = 3 element vector. 24-bit RGB representation of the average
%               color of sand.
%
% Keith Ma, July 2015

% DEBUG: Run as a script with predefined inputs
image_file = 'test_data/test_a_00.jpg';

%% read in data, check for sane inputs

