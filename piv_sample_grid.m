function [r_sample, c_sample] = ...
             piv_sample_grid(sample_spc, x_world, y_world)
% function [r_sample, c_sample, x_sample_vector, y_sample_vector] = ...
%              piv_sample_grid(sample_spc, x_world, y_world)
%
% Create sample grid which is the largest centered grid that fits in the domain
%
% Arguments:
%   sample_spc = Scalar, integer, grid spacing in pixels
%   x_world, y_world = Vectors, double, world coordinates for the input images.
%       Vector sizes match the input images, such that size(image, 1) ==
%       length(y_world) and size(image, 2) == length(x_world)
%   r_sample, c_sample = 2D matrix, integer, coordinate matrices, as constructed
%       by meshgrid, for the sample grid in the y (a.k.a. row) and x (a.k.a.
%       column) directions, in pixels
% %

remainder = mod((length(y_world)-1), sample_spc);
r_sample_vector = (1 + remainder/2):sample_spc:length(y_world);

remainder = mod((length(x_world)-1), sample_spc);
c_sample_vector = (1 + remainder/2):sample_spc:length(x_world);

[c_sample, r_sample] = meshgrid(c_sample_vector, r_sample_vector);

% % <DEBUG>: Check grid edges to confirm it is centered in the domain
% fprintf('Left edge [pixel]: %f, Right edge [pixel]: %f\n', ...
%     c_sample_vector(1) - 1, length(x_world) - c_sample_vector(end));
% fprintf('Left edge [world]: %f, Right edge [world]: %f\n', ...
%     x_sample_vector(1) - x_world(1), x_world(end) - x_sample_vector(end));
% fprintf('Bottom edge [pixel]: %f, Top edge [pixel]: %f\n', ...
%     r_sample_vector(1) - 1, length(y_world) - r_sample_vector(end));
% fprintf('Bottom edge [world]: %f, Top edge [world]: %f\n', ...
%     y_sample_vector(1) - y_world(1), y_world(end) - y_sample_vector(end));
% % </DEBUG>