function [xx, yy, uu_piv, vv_piv, uu_exact, vv_exact] = test_piv_v2(varargin)
% function [] = test_piv(varargin)
%
% Test PIV performance using case using synthetic image pair
%
% NOTE: This is a "high-level" test that covers all of the PIV subsystems.
%
% Optional Arguments ('Name', Value):
%   'image_file': String, name of pre-processed image file, as produced by
%       prep_series(), to read for raw image, default = ./test/default_image.nc
%   'image_index': Integer, (1-based) index of image in image file to use for
%       raw image, default = 1
%   'image_pos': 4-element position vector indicating the limits of the image to
%       extract and deform. Must contain only sand (all within the ROI),
%       in meters, default = [-0.12, 0.005, 0.092, 0.07]
%   'u1': 2-element vector specifying first end-member translation in pixels,
%       default = 20*[cosd(45), sind(45)]
%   'u2': 2-element vector specifying second end-member translation in pixels,
%       default = 10*[cosd(45), sind(45)]
%   'theta': Scalar, orientation of shear band specified as
%       counter-clockwise angle to the positive x-axis, in degrees, limited to
%       range 0 - 90, default = 45
%   'shear_width': Scalar, width of shear band, as fraction of image width,
%       default = 0.25   
%   'pad_width': Scalar, width of edge padding to add to image (to accomodate
%       edge displacements) as a fraction of image size, default = 0.1
%   'samp_len': piv() parameter, default [60, 30]
%   'samp_spc': piv() parameter, default 15
%   'intr_len': piv() parameter, default [90, 40]
%   'num_pass': piv() parameter, default [1, 2]
%   'valid_radius: piv() parameter, default 45
%   'valid_max': piv() parameter, default 2
%   'valid_eps': piv() parameter, default 0.1
%   'min_frac_data': piv() parameter, default 0.5
%   'min_frac_overlap': piv() parameter, default 0.25
%   'verbose': Scalar logical, set true to enable verbose reporting for all
%       components of the analysis
% %

