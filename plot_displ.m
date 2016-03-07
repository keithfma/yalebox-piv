function opt = plot_displ(piv_file, index, varargin)
% function opt = plot_displ(piv_file, index, varargin)
%
% Plot displacement magnitude and direction. Generates 3 subplots:
%
%   - x-direction displacement 
%   - y-direction displacement 
%   - total displacement magnitude and direction
% 
% Arguments:
%
%   piv_file = String, name of PIV results netCDF file, as produced by
%       piv_series()
%
%   index = Scalar, integer, index of time-step in piv_file to plot
% 
%   opt = Struct, options used for this plotting run. Primarily useful for
%       plotting series.
%   
% Parameters (Name-Value pairs):
%
%   'coord_units' = String, name of units for coordinate axes, coordinate values
%       will be rescaled accordingly, and labels will reflect these units. Valid
%       options are: 'm', 'cm'. Default = 'cm'
%
%   'displ_units' = String, name of units for displacement vector components and
%       magnitude, values will be rescaled accordingly, and labels will reflect
%       these units. Valid options are: 'm', 'mm', '1'. Default is 'mm'.
%
%   'norm_bbox' = Vector, length==4, bounding box for the data region used in
%       computing the displacement unit normalization factor. This is only used
%       if displ_units == '1'. If empty, the normalization function will prompt
%       the user to select a box interactively. Default = [].
%
%   'xlim', 'ylim', 'ulim', 'vlim', 'mlim' = Vector, length==2, [minimum,
%       maximum] values for the x-axis, y-axis, u-displacement, v-displacement,
%       and displacement magnitude. Values beyond the range of the data will be
%       truncated to the data limits. Thus, to span the data, one could use
%       [-inf, inf].
%
%   'qsize' = Vector, length==2, size of the grid for vector direction overlay
%       (quiver) in [rows, cols]
%
%   'qbnd' = Scalar, boundary margin for vector direction overlay, as fraction
%       of the axis ranges.
%
%   'qscale' = Scalar, range [0,1], length of vector lines in vector direction
%       overlay 
%
%   'font_size_title', 'font_size_tick', 'font_size_axis' = Scalar,
%       integer, font size for titles, tick labels, and axis labels,
%       repectively. Defaults are 14, 12, 12.
% %

%% parse input arguments

% constant parameters
spline_tension = 0.9;

% positional arguments
validateattributes(piv_file, {'char'}, {'vector'});
assert(exist(piv_file, 'file')==2);
validateattributes(index, {'numeric'}, {'scalar', 'integer'});

% parameter name-value pairs
ip = inputParser();

ip.addParameter('coord_units', 'cm', ...
    @(x) ismember(x, {'m', 'cm'})); 
ip.addParameter('displ_units', 'mm/step', ...
    @(x) ismember(x, {'m/step', 'mm/step', '1'}));
ip.addParameter('norm_bbox', [], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 4}));
ip.addParameter('xlim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('ylim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('ulim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('vlim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('mlim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('qsize', [20, 10], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2, 'integer'}));
ip.addParameter('qbnd', 0.05, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 0.5}));
ip.addParameter('qscale', 0.2, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
ip.addParameter('font_size_title', 14, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
ip.addParameter('font_size_tick', 12, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
ip.addParameter('font_size_axis', 12, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));

ip.parse(varargin{:});
opt = ip.Results;

%% prepare data

[step, xx, yy, uu, vv, mm] = util_read_piv_step(piv_file, index);

[xx, yy, uu, vv, mm, opt.norm_bbox] = ...
    util_convert_displ_units(xx, yy, uu, vv, mm, opt.coord_units, ...
        opt.displ_units, opt.norm_bbox);

% truncate axis limits to data range
opt.xlim(1) = max(opt.xlim(1), min(xx)); 
opt.ylim(1) = max(opt.ylim(1), min(yy));
opt.ulim(1) = max(opt.ulim(1), min(uu(:)));
opt.vlim(1) = max(opt.vlim(1), min(vv(:))); 
opt.mlim(1) = max(opt.mlim(1), min(mm(:))); 
opt.xlim(2) = min(opt.xlim(2), max(xx));
opt.ylim(2) = min(opt.ylim(2), max(yy)); 
opt.ulim(2) = min(opt.ulim(2), max(uu(:)));
opt.vlim(2) = min(opt.vlim(2), max(vv(:)));
opt.mlim(2) = min(opt.mlim(2), max(mm(:)));

% quiver plot data: interpolate to low-res checkerboard grid, normalize length
%...generate staggered regular grid
qx0 = opt.xlim(1) + range(opt.xlim)*opt.qbnd;
qy0 = opt.ylim(1) + range(opt.ylim)*opt.qbnd;
qx1 = opt.xlim(2) - range(opt.xlim)*opt.qbnd;
qy1 = opt.ylim(2) - range(opt.ylim)*opt.qbnd;
qx = linspace(qx0, qx1, opt.qsize(2));
qy = linspace(qy0, qy1, opt.qsize(1));
[qx, qy] = meshgrid(qx, qy);
chkr = bsxfun(@xor, mod(1:opt.qsize(2), 2), mod(1:opt.qsize(1), 2)');
qx = qx(chkr);
qy = qy(chkr);
%...crop grid to the roi
[xgrid, ygrid] = meshgrid(xx, yy);
roi = ~isnan(uu);
roiq = interp2(xgrid, ygrid, roi, qx, qy, 'nearest');
qx = qx(roiq);
qy = qy(roiq);
%...interpolate vectors
qu = spline2d(qx, qy, xgrid(roi), ygrid(roi), uu(roi), spline_tension);
qv = spline2d(qx, qy, xgrid(roi), ygrid(roi), vv(roi), spline_tension);
%...normalize vector length
mq = sqrt(qu.^2+qv.^2);
qu = qu./mq;
qv = qv./mq;

%% plot

[ax_top, ax_mid, ax_bot] = create_figure();

% plot displacement magnitude and direction
axes(ax_top);
imagesc(xx, yy, mm, 'AlphaData', ~isnan(mm)); 
plot_direction(qx, qy, qu, qv, opt.qscale);
format_axes(gca, opt, 0, 1, opt.mlim, ...
    sprintf('displacement vector magnitude, step = %.1f', step));

% plot x-direction displacement magnitude
axes(ax_mid);
imagesc(xx, yy, uu, 'AlphaData', ~isnan(uu)); 
plot_direction(qx, qy, qu, qv, opt.qscale);
format_axes(gca, opt, 0, 1, opt.ulim, ...
    sprintf('x-direction displacement component, step = %.1f', step));

% plot y-direction displacement magnitude
axes(ax_bot);
imagesc(xx, yy, vv, 'AlphaData', ~isnan(vv)); 
plot_direction(qx, qy, qu, qv, opt.qscale);
format_axes(gca, opt, 1, 1, opt.vlim, ...
    sprintf('y-direction displacement component, step = %.1f', step));

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

function [] = plot_direction(qx, qy, qu, qv, qscale)
% overlay a simple vector plot

hold on
quiver(qx, qy, qu, qv, qscale, 'Color', 'w');
hold off


function [] = format_axes(ax, opt, xshow, yshow, clim, title_str)
% apply consistent formatting to axes

axis equal
grid on
set(ax, 'XLim', opt.xlim, ...
        'YLim', opt.ylim, ...
        'YDir', 'normal', ...
        'Color', 0.9*[1, 1, 1], ...
        'FontSize', opt.font_size_tick);
if xshow
    xlabel(sprintf('x [%s]', opt.coord_units), 'FontSize', opt.font_size_axis);
else
    set(ax, 'XTickLabel', '');
end
if yshow
    ylabel(sprintf('y [%s]', opt.coord_units), 'FontSize', opt.font_size_axis);
else
    set(ax, 'YTickLabel', '');
end   
caxis(clim);
h = colorbar;
ylabel(h, sprintf('displacement [%s]', opt.displ_units), ...
    'FontSize', opt.font_size_axis);
title(title_str, 'FontSize', opt.font_size_title);