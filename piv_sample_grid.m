function [r_sample, c_sample] = piv_sample_grid(sample_spc, img_nr, img_nc)
%function [r_sample, c_sample] = piv_sample_grid(sample_spc, img_nr, img_nc)
%
% Create sample grid which is the largest centered grid that fits in the domain
%
% Arguments:
%   sample_spc = Scalar, integer, grid spacing in pixels
%   img_nr, img_nc = Scalar, dimensions (number of rows and columns) of the
%       input images.
%   r_sample, c_sample = 2D matrix, integer, coordinate matrices, as constructed
%       by meshgrid, for the sample grid in the y (a.k.a. row) and x (a.k.a.
%       column) directions, in pixels
% %

remainder = mod((img_nr - 1), sample_spc);
r_sample_vector = (1 + remainder/2):sample_spc:img_nr;

remainder = mod((img_nc - 1), sample_spc);
c_sample_vector = (1 + remainder/2):sample_spc:img_nc;

[c_sample, r_sample] = meshgrid(c_sample_vector, r_sample_vector);

% % <DEBUG>: Check grid edges to confirm it is centered in the domain
% fprintf('Left edge [pixel]: %f, Right edge [pixel]: %f\n', ...
%     c_sample_vector(1) - 1, img_nc - c_sample_vector(end));
% fprintf('Left edge [world]: %f, Right edge [world]: %f\n', ...
%     x_sample_vector(1) - x_world(1), x_world(end) - x_sample_vector(end));
% fprintf('Bottom edge [pixel]: %f, Top edge [pixel]: %f\n', ...
%     r_sample_vector(1) - 1, img_nr - r_sample_vector(end));
% fprintf('Bottom edge [world]: %f, Top edge [world]: %f\n', ...
%     y_sample_vector(1) - y_world(1), y_world(end) - y_sample_vector(end));
% % </DEBUG>