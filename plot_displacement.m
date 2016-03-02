function [] = plot_displacement(piv_file, index, xlimits, ylimits, bbox)
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
% xlimits, ylimits = Vector, length==2, [minimum, maximum] values for the x- and
%   y-axes, values larger than the range of the input data will be truncated.
%   Thus, to span the data, one could use [-inf, inf].
%
% bbox = (optional) Vector, length==4, bounding box for the data region to be
%   used to compute displacement normalization. Provided in world coordinates,
%   formatted as [left, bottom, width, height]. *If empty, data are not
%   normalized*
% %

% debug: hard-code input arguments
piv_file = '~/Documents/dissertation/yalebox-exp-fault/data/fault_ss_01/piv/fault_ss_01_sidef.displ.nc';
index = 100;
xlimits = [-inf, inf];
ylimits = [-inf, 0.1];
bbox = [];

% normalize?
normalize = true;
if isempty(bbox)
    normalize = false;
end

% check for sane inputs
validateattributes(piv_file, {'char'}, {'vector'});
validateattributes(index, {'numeric'}, {'scalar', 'positive', 'integer'});
validateattributes(xlimits, {'numeric'}, {'vector', 'numel', 2, 'increasing'});
validateattributes(ylimits, {'numeric'}, {'vector', 'numel', 2, 'increasing'});
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

% truncate axis limits
xlimits(1) = max(xlimits(1), min(xx));
xlimits(2) = min(xlimits(2), max(xx));
ylimits(1) = max(ylimits(1), min(yy));
ylimits(2) = min(ylimits(2), max(yy));

% convert units
if ~normalize
    uu = uu*1000; % m -> mm
    vv = vv*1000; % m -> mm
    mm = mm*1000; % m -> mm
    color_units = '[mm/step]';
end

% create figures and axes
[ax_top, ax_mid, ax_bot] = create_figure();

% plot displacement magnitude and direction
axes(ax_top)
imagesc(xx, yy, mm, 'AlphaData', ~isnan(mm)); 
format_axes(gca, xlimits, ylimits, color_units);

% plot x-direction displacement magnitude
axes(ax_mid);
imagesc(xx, yy, uu, 'AlphaData', ~isnan(uu)); 
format_axes(gca, xlimits, ylimits, color_units);

% plot y-direction displacement magnitude
axes(ax_bot);
imagesc(xx, yy, vv, 'AlphaData', ~isnan(vv)); 
format_axes(gca, xlimits, ylimits, color_units);

keyboard

end

function [top, mid, bot] = create_figure()
% setup figure with 3 x 1 subplots with tight spacing

% constant parameters
figure_position = [0, 0, 1, 1]; 
spc_left = 0.1;
spc_right = 0.1;
spc_top = 0.1;
spc_mid = 0.1; 
spc_bot = 0.1;

% derived parameters
height_ax = (1-spc_top-2*spc_mid-spc_bot)/3;  
width_ax = 1-spc_left-spc_right;
bot_position = [spc_left, spc_bot,                       width_ax, height_ax];
mid_position = [spc_left, spc_bot+height_ax+spc_mid,     width_ax, height_ax];
top_position = [spc_left, spc_bot+2*height_ax+2*spc_mid, width_ax, height_ax];

figure
set(gcf, 'Units', 'Normalized', 'Position', figure_position);
bot = axes('Position', bot_position);
mid = axes('Position', mid_position);
top = axes('Position', top_position);

end

function [] = format_axes(ax, xlimits, ylimits, color_units)

% constant parameters
axes_color = 0.9*[1, 1, 1];

% format
axis equal
h = colorbar;
ylabel(h, color_units);
set(ax, 'XLim', xlimits, ...
        'YLim', ylimits, ...
        'YDir', 'normal', ...
        'Color', axes_color);
    
end
