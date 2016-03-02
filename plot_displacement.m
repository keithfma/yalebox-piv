function [] = plot_displacement(piv_file, index, bbox)
% function [] = plot_displacement(xx, yy, uu, vv, bbox)
%
% Plot normalized displacement magnitude and direction. Function generates 3 separate plots:
% - x-direction displacement 
% - y-direction displacement 
% - total displacement magnitude and direction
% 
% Arguments:
%
% piv_file = String, file name of netCDF file containing output from the
%   piv_series() analysis
%
% index = Scalar, index of timestep in piv_file to plot, uses MATLAB-style
%   1-based indices.
%
% bbox = (optional) Vector, length==4, bounding box for the data region to be
%   used to compute displacement normalization. Provided in world coordinates,
%   formatted as [left, bottom, width, height]. *If empty, data are not
%   normalized*
% %

% debug: hard-code input arguments
piv_file = '~/Documents/dissertation/yalebox-exp-fault/data/fault_ss_01/piv/fault_ss_01_sidef.displ.nc';
index = 100;
bbox = [];

% local constants
font_size = 14;

% normalize?
normalize = true;
if isempty(bbox)
    normalize = false;
end

% check for sane inputs
validateattributes(piv_file, {'char'}, {'vector'});
validateattributes(index, {'numeric'}, {'scalar', 'positive', 'integer'});
if normalize
    validateattributes(bbox, {'numeric'}, {'vector', 'real', 'numel', 4});
end

% read data from netCDF
xx = ncread(piv_file, 'x'); nx = numel(xx);
yy = ncread(piv_file, 'y'); ny = numel(yy);
step = ncread(piv_file, 'step', index, 1);
uu = squeeze(ncread(piv_file, 'u', [1, 1, index], [nx, ny, 1]))';
vv = squeeze(ncread(piv_file, 'v', [1, 1, index], [nx, ny, 1]))';