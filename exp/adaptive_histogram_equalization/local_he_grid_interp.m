function [eql] = local_he_grid_interp(im, ignore)
%
% Local (adaptive) histogram equalization using brute-force approach.
% Equalization tranform is computed for each pixel based on a neighboring
% pixels.
%
% Arguments:
%
% im = 2D matrix, double. Continuous-valued image matrix.
%
% ignore = Scalar, double. Values to ignore in equalization routine (masked
%   pixels).
% 
% eql = 2D matrix, double. Equalized image matrix, with a uniform distribution
%   in the range [0,1]
%
% %

% check inputs
validateattributes(im, {'double'}, {'2d', 'real'}, mfilename, 'im');
validateattributes(ignore, {'double'}, {'scalar', 'real'}, mfilename, 'ignore');

% dummy output
eql = im;