function [eql] = local_he_brute_force(im, ignore, nwin)
% function [eql] = local_he_brute_force(im, ignore, nwin)
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
% nwin = Scalar, integer, odd. Side length (in pixels) for the local
%   neighborhood used to compute the transform for each pixel.
%
% eql = 2D matrix, double. Equalized image matrix, with a uniform distribution
%   in the range [0,1]
%
% %

% check inputs
validateattributes(im, {'double'}, {'2d', 'real'}, mfilename, 'im');
validateattributes(ignore, {'double'}, {'scalar', 'real'}, mfilename, 'ignore');
validateattributes(nwin, {'numeric'}, {'integer', 'positive', 'odd'}, mfilename, 'nwin');

% debug: dummy output
eql = zeros(size(im));
% } debug