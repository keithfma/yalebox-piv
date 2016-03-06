function plot_displ(step, xx, yy, uu, vv, mm, pp)
% function plot_displ(step, xx, yy, uu, vv, mm, pp)
%
% Plot displacement magnitude and direction. Generates 3 subplots:
%   - x-direction displacement 
%   - y-direction displacement 
%   - total displacement magnitude and direction
% 
% Arguments:
%
% xx, yy = 
%
% uu, vv, mm = 
%
% pp = Struct, contains all parameters needed for plot formating, etc. 
%   Members:
%
%       .xunit, .yunit, .uunit, .vunit, .munit = String, unit names for the
%           x-axis, y-axis, u-displacement, v-displacement, and displacement
%           magnitude. Used to generate labels.
%
%       .xlim, .ylim, .ulim, .vlim, .mlim = Vector, length==2, [minimum,
%           maximum] values for the x-axis, y-axis, u-displacement,
%           v-displacement, and displacement magnitude. Values beyond the range
%           of the data will be truncated to the data limits. Thus, to span the
%           data, one could use [-inf, inf].
%
%       .qsize = Vector, length==2, size of the grid for vector direction
%           overlay (quiver) in [rows, cols].
%
%       .qbnd = Scalar, boundary margin for vector direction overlay, as
%           fraction of the axis ranges.
%
%       .qscale = Scalar, range [0,1], length of vector lines in vector
%           direction overlay
% %

% set defaults
if nargin < 7  || isempty(pp)
    pp = default_plot_displ;
end

% local parameters
tension = 0.9;
            
% truncate limits to data
pp.xlim(1) = max(pp.xlim(1), min(xx)); 
pp.ylim(1) = max(pp.ylim(1), min(yy));
pp.ulim(1) = max(pp.ulim(1), min(uu(:)));
pp.vlim(1) = max(pp.vlim(1), min(vv(:))); 
pp.mlim(1) = max(pp.mlim(1), min(mm(:))); 
pp.xlim(2) = min(pp.xlim(2), max(xx));
pp.ylim(2) = min(pp.ylim(2), max(yy)); 
pp.ulim(2) = min(pp.ulim(2), max(uu(:)));
pp.vlim(2) = min(pp.vlim(2), max(vv(:)));
pp.mlim(2) = min(pp.mlim(2), max(mm(:)));

% interpolate vectors to low-res checkerboard grid, normalize length for quiver plots
%...generate staggered regular grid
qx = linspace(xlim(1)+range(xlim)*pp.qbnd, xlim(2)-range(xlim)*pp.qbnd, pp.qsize(2));
qy = linspace(ylim(1)+range(ylim)*pp.qbnd, ylim(2)-range(ylim)*pp.qbnd, pp.qsize(1));
[qx, qy] = meshgrid(qx, qy);
chkr = bsxfun(@xor, mod(1:pp.qsize(2), 2), mod(1:pp.qsize(1), 2)');
qx = qx(chkr);
qy = qy(chkr);
%...crop grid to the roi
[xgrid, ygrid] = meshgrid(xx, yy);
roi = ~isnan(uu);
roiq = interp2(xgrid, ygrid, roi, qx, qy, 'nearest');
qx = qx(roiq);
qy = qy(roiq);
%...interpolate vectors
qu = spline2d(qx, qy, xgrid(roi), ygrid(roi), uu(roi), tension);
qv = spline2d(qx, qy, xgrid(roi), ygrid(roi), vv(roi), tension);
%...normalize vector length
mq = sqrt(qu.^2+qv.^2);
qu = qu./mq;
qv = qv./mq;

% create figures and axes
[ax_top, ax_mid, ax_bot] = create_figure();

% plot displacement magnitude and direction
axes(ax_top)
imagesc(xx, yy, mm, 'AlphaData', ~isnan(mm)); 
plot_direction(qx, qy, qu, qv, pp.qscale);
format_axes(gca, pp.xlim, pp.xunits, 0, pp.ylim, pp.yunits, 1, pp.mlim, pp.munits, ...
    sprintf('displacement vector magnitude, step = %.1f', step));

% plot x-direction displacement magnitude
axes(ax_mid);
imagesc(xx, yy, uu, 'AlphaData', ~isnan(uu)); 
plot_direction(qx, qy, qu, qv, pp.qscale);
format_axes(gca, pp.xlim, pp.xunits, 0, pp.ylim, pp.yunits, 1, pp.ulim, pp.uunits, ...
    sprintf('x-direction displacement component, step = %.1f', step));

% plot y-direction displacement magnitude
axes(ax_bot);
imagesc(xx, yy, vv, 'AlphaData', ~isnan(vv)); 
plot_direction(qx, qy, qu, qv, qscale);
format_axes(gca, pp.xlim, pp.xunits, 1, pp.ylim, pp.yunits, 1, pp.vlim, pp.vunits, ...
    sprintf('y-direction displacement component, step = %.1f', step));

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

function [] = plot_direction(qx, qy, qu, qv, qscale)
% overlay a simple vector plot

hold on
quiver(qx, qy, qu, qv, qscale, 'Color', 'w');
hold off

end

function [] = format_axes(ax, xlim, xunits, xshow, ylim, yunits, yshow, ...
                climits, cunits, title_str)
% apply consistent formatting to axes

axis equal
grid on
set(ax, 'XLim', xlim, ...
        'YLim', ylim, ...
        'YDir', 'normal', ...
        'Color', 0.9*[1, 1, 1]);
if xshow
    xlabel(['x ', xunits]); 
else
    set(ax, 'XTickLabel', '');
end
if yshow
    ylabel(['y ', yunits]); 
else
    set(ax, 'YTickLabel', '');
end   
caxis(climits);
h = colorbar;
ylabel(h, cunits);
title(title_str);
    
end