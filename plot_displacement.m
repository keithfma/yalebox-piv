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

% local parameters 
font_size = 14;
figure_position = [0, 0, 1, 1]; % normalized, [L, B, W, H]
spc_left = 0.1;
spc_right = 0.1;
spc_top = 0.1;
spc_mid = 0.1; 
spc_bot = 0.1;

% local derived parameters
height_ax = (1-spc_top-2*spc_mid-spc_bot)/3;  
width_ax = 1-spc_left-spc_right;
ax_bot_position = [spc_left, spc_bot,                       width_ax, height_ax];
ax_mid_position = [spc_left, spc_bot+height_ax+spc_mid,     width_ax, height_ax];
ax_top_position = [spc_left, spc_bot+2*height_ax+2*spc_mid, width_ax, height_ax];

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
mm = sqrt(uu.*uu+vv.*vv);

% setup figure with 3 x 1 subplots with tight spacing
figure
set(gcf, 'Units', 'Normalized', 'Position', figure_position);
ax_bot = axes('Position', ax_bot_position);
ax_mid = axes('Position', ax_mid_position);
ax_top = axes('Position', ax_top_position);

% plot displacement magnitude and direction
axes(ax_top)
imagesc(xx, yy, mm, 'AlphaData', ~isnan(mm)); 
set(gca, 'YDir', 'normal');

% plot x-direction displacement magnitude
axes(ax_mid);
imagesc(xx, yy, uu, 'AlphaData', ~isnan(uu)); 
set(gca, 'YDir', 'normal');

% plot y-direction displacement magnitude
axes(ax_bot);
imagesc(xx, yy, vv, 'AlphaData', ~isnan(vv)); 
set(gca, 'YDir', 'normal');

keyboard