function [eql] = local_he_grid_interp(im, len, spc, ignore)
%
% Local (adaptive) histogram equalization using gridded interpolation approach.
% Equalization tranforms are computed for on a low resolution grid, and pixels
% values are assigned by bilinear interpolation of the four neighboring
% transforms. Edges are treated by padding the transform grid .
%
% Arguments:
%
% im = 2D matrix, double. Continuous-valued image matrix.
%
% ignore = Scalar, double. Values to ignore in equalization routine (masked
%   pixels).
%
% len = Scalar, integer, odd. Side length (in pixels) for the local
%   neighborhood used to compute the transform for each pixel.
% 
% eql = 2D matrix, double. Equalized image matrix, with a uniform distribution
%   in the range [0,1]
%
% %

% check inputs
validateattributes(im, {'double'}, {'2d', 'real'}, mfilename, 'im');
validateattributes(ignore, {'double'}, {'scalar', 'real'}, mfilename, 'ignore');

% ABANDONED THIS ATTEMPT FOR NOW
error('local_he_grid_interp not implemented');

% dummy output
eql = im;