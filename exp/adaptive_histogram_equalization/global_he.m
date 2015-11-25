function [eql] = global_he(im)
%
% Remap intensities of the input image so that the intensity PDF is (nearly)
% uniform.
%
% im = 2D matrix, double. Continuous-valued image matrix.
%
% eql = 2D matrix, double. Equalized image matrix, with a uniform distribution
% in the range [0,1]

validateattributes(im, {'double'}, {'2d', 'real'}, mfilename, 'im');

% debug: set output variable {
eql = im; 
% } debug