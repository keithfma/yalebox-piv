function [] = plot_displacement(xx, yy, uu, vv, bbox)
% function [] = plot_displacement(xx, yy, uu, vv, bbox)
%
% Plot normalized displacement magnitude and direction. Function generates 3 separate plots:
% - x-direction displacement 
% - y-direction displacement 
% - total displacement magnitude and direction
% 
% Arguments:
%
% xx, yy = Vectors, x- and y-direction coordinate vectors in [meters] with the
%   origin at the s-point, corresponding to uu, vv
%
% uu, vv = 2D matrix, x- and y- direction displacement vector components in
%   [meters/step]
%
% bbox = Vector, length==4, bounding box for the data region to be used to
%       compute displacement normalization. Provided in world coordinates, e.g.
%       [meters] with the origin at the s-point, formatted as 
%       [left, bottom, width, height]
% %

% local constants