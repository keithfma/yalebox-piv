function [bbox] = plot_displacement_norm(piv_file, index, bbox, xlim, ylim, ulim, vlim, mlim, ...
                  qsize, qbnd, qscale)
%
% Plot displacement magnitude and direction. Magnitudes are normalised to
% facilitate comparison between steps of unequal size. This inconsistency is a
% shortcoming of the experimental apparatus. The normalization factor is a
% user-selected characteristic velocity, computed as the absolute value of the
% median of the velocity mangitude in a user-defined window. For most sandbox
% experiments, the logical choice is the incoming section.
% 
% Function generates 3 separate plots:
%
%   - x-direction displacement 
%   - y-direction displacement 
%   - total displacement magnitude and direction
%
% Basic options to customize the plot are included as optional arguments, others
% are hard-coded.
% 
% Arguments:
%
% piv_file = String, file name of netCDF file containing output from the
%   piv_series() analysis
%
% index = Scalar, index of timestep in piv_file to plot, uses MATLAB-style
%   1-based indices.
%
% bbox = (Optional) Vector, length==4, bounding box for the data region used to
%   normalize velocity mangitudes. The expected format is [x0, y0, xwidth,
%   yheight] in world coordinates, where x0, y0 define the lower left corner of
%   the box, and xwidth, yheight define its size. If the argument is left empty,
%   the box must be selected interactively. Due to this interactive option, the
%   value of bbox is returned for reuse.
%
% xlim, ylim, ulim, vlim, mlim = (Optional) Vector, length==2, [minimum,
%   maximum] values for the x-axis, y-axis, u-displacement, v-displacement, and
%   displacement magnitude. Values beyond the range of the data will be
%   truncated to the data limits. Thus, to span the data, one could use [-inf, inf].
%
% qsize = (Optional) Vector, length==2, size of the grid for vector direction
%   overlay (quiver) in [rows, cols]. 
%
% qbnd = (Optional) Scalar, boundary margin for vector direction overlay, as
%   fraction of the axis ranges.
%
% qscale = (Optional) Scalar, range [0,1], length of vector lines in vector
%   direction overlay 
% %

% set defaults
if nargin < 4  || isempty(xlim);   xlim = [-inf, inf]; end
if nargin < 5  || isempty(ylim);   ylim = [-inf, inf]; end
if nargin < 6  || isempty(ulim);   ulim = [-inf, inf]; end
if nargin < 7  || isempty(vlim);   vlim = [-inf, inf]; end
if nargin < 8  || isempty(mlim);   mlim = [-inf, inf]; end
if nargin < 9  || isempty(qsize);  qsize = [20, 10];   end
if nargin < 10 || isempty(qbnd);   qbnd = 0.05;        end
if nargin < 11 || isempty(qscale); qscale = 0.1;       end

% check for sane inputs
validateattributes(piv_file, {'char'}, {'vector'});
assert(exist(piv_file, 'file') == 2);
validateattributes(index, {'numeric'}, {'scalar', 'integer'});
if ~isempty(bbox)
    validateattributes(bbox, {'numeric'}, {'vector', 'numel', 4});
end
validateattributes(xlim, {'numeric'}, {'vector', 'numel', 2});
validateattributes(ylim, {'numeric'}, {'vector', 'numel', 2});
validateattributes(ulim, {'numeric'}, {'vector', 'numel', 2});
validateattributes(vlim, {'numeric'}, {'vector', 'numel', 2});
validateattributes(mlim, {'numeric'}, {'vector', 'numel', 2});
validateattributes(qsize, {'numeric'}, {'vector', 'numel', 2, 'integer'});
validateattributes(qbnd, {'numeric'}, {'scalar', '>=', 0, '<=', 1});
validateattributes(qscale, {'numeric'}, {'scalar', '>=', 0, '<=', 1});

% local parameters
tension = 0.9;

% read data from netCDF
xx = ncread(piv_file, 'x'); 
xx = double(xx);
nx = numel(xx);

yy = ncread(piv_file, 'y'); 
yy = double(yy);
ny = numel(yy);

step = ncread(piv_file, 'step', index, 1);
step = double(step);

uu = ncread(piv_file, 'u', [1, 1, index], [nx, ny, 1])';
uu = double(squeeze(uu));

vv = ncread(piv_file, 'v', [1, 1, index], [nx, ny, 1])';
vv = double(squeeze(vv));

mm = sqrt(uu.*uu+vv.*vv);

% interactively select bounding box for normalization, if needed
if isempty(bbox);    
    figure;
    imagesc(xx, yy, mm, 'AlphaData', ~isnan(mm));
    set(gca, 'YDir', 'Normal');
    title('Select bounding box for normalization using mouse');
    bbox = getrect();
    close(gcf);
end

% get normalization factor: median displacement magnitude within bounding box
x_in_bbox = (xx >= bbox(1)) & (xx <= (bbox(1)+bbox(3)));
y_in_bbox = (yy >= bbox(2)) & (yy <= (bbox(2)+bbox(4)));
in_bbox = logical(bsxfun(@times, x_in_bbox(:)', y_in_bbox(:)));
norm = median(mm(in_bbox));
            
% truncate axis limits
xlim(1) = max(xlim(1), min(xx)); xlim(2) = min(xlim(2), max(xx));
ylim(1) = max(ylim(1), min(yy)); ylim(2) = min(ylim(2), max(yy));
ulim(1) = max(ulim(1), min(uu(:))); ulim(2) = min(ulim(2), max(uu(:)));
vlim(1) = max(vlim(1), min(vv(:))); vlim(2) = min(vlim(2), max(vv(:)));
mlim(1) = max(mlim(1), min(mm(:))); mlim(2) = min(mlim(2), max(mm(:)));

% convert units
xunits = '[cm]';
xx = xx*100; % m -> cm
xlim = xlim*100;

yunits = '[cm]';
yy = yy*100; % m -> cm
ylim = ylim*100;

cunits = '[1]';
uu = uu/norm; % m -> mm
vv = vv/norm;
mm = mm/norm;
ulim = ulim/norm;
vlim = vlim/norm;
mlim = mlim/norm;

% interpolate vectors to low-res checkerboard grid, normalize length for quiver plots
%...generate staggered regular grid
qx = linspace(xlim(1)+range(xlim)*qbnd, xlim(2)-range(xlim)*qbnd, qsize(2));
qy = linspace(ylim(1)+range(ylim)*qbnd, ylim(2)-range(ylim)*qbnd, qsize(1));
[qx, qy] = meshgrid(qx, qy);
chkr = bsxfun(@xor, mod(1:qsize(2), 2), mod(1:qsize(1), 2)');
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
plot_direction(qx, qy, qu, qv, qscale);
format_axes(gca, xlim, xunits, 0, ylim, yunits, 1, mlim, cunits, ...
    sprintf('displacement vector magnitude, step = %.1f', step));

% plot x-direction displacement magnitude
axes(ax_mid);
imagesc(xx, yy, uu, 'AlphaData', ~isnan(uu)); 
plot_direction(qx, qy, qu, qv, qscale);
format_axes(gca, xlim, xunits, 0, ylim, yunits, 1, ulim, cunits, ...
    sprintf('x-direction displacement component, step = %.1f', step));

% plot y-direction displacement magnitude
axes(ax_bot);
imagesc(xx, yy, vv, 'AlphaData', ~isnan(vv)); 
plot_direction(qx, qy, qu, qv, qscale);
format_axes(gca, xlim, xunits, 1, ylim, yunits, 1, vlim, cunits, ...
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