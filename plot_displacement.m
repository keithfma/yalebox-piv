function [] = plot_displacement(piv_file, index, xlim, ylim, ulim, vlim, mlim)
%
% Plot displacement magnitude and direction. Function generates 3 separate plots:
%   - x-direction displacement 
%   - y-direction displacement 
%   - total displacement magnitude and direction
% Converts units to [cm] for axes and to [mm/step] for displacements
% 
% Arguments:
%
% piv_file = String, file name of netCDF file containing output from the
%   piv_series() analysis
%
% index = Scalar, index of timestep in piv_file to plot, uses MATLAB-style
%   1-based indices.
%
% xlim, ylim, ulim, vlim, mlim = Vector, length==2, [minimum, maximum] values for the x-axis, y-axis, u-displacement, v-displacement, and
%   displacement magnitude. Values beyond the range of the data will be
%   truncated to the data limits. Thus, to span the data, one could use [-inf, inf].
% %

% set defaults
if nargin < 3 || isempty(xlim); xlim = [-inf, inf]; end
if nargin < 4 || isempty(ylim); ylim = [-inf, inf]; end
if nargin < 5 || isempty(ulim); ulim = [-inf, inf]; end
if nargin < 6 || isempty(vlim); vlim = [-inf, inf]; end
if nargin < 7 || isempty(mlim); mlim = [-inf, inf]; end

% local parameters
tension = 0.9;

% debug: hard-code input arguments
piv_file = '~/Documents/dissertation/yalebox-exp-fault/data/fault_ss_01/piv/fault_ss_01_sidef.displ.nc';
index = 250;
xlim = [-inf, inf];
ylim = [-inf, 0.12];
ulim = [-inf, inf];
vlim = [-inf, inf];
mlim = [-inf, inf];
% end debug

% check for sane inputs
validateattributes(piv_file, {'char'}, {'vector'});
validateattributes(index, {'numeric'}, {'scalar', 'positive', 'integer'});
validateattributes(xlim, {'numeric'}, {'vector', 'numel', 2, 'increasing'});
validateattributes(ylim, {'numeric'}, {'vector', 'numel', 2, 'increasing'});
validateattributes(ulim, {'numeric'}, {'vector', 'numel', 2, 'increasing'});
validateattributes(vlim, {'numeric'}, {'vector', 'numel', 2, 'increasing'});
validateattributes(mlim, {'numeric'}, {'vector', 'numel', 2, 'increasing'});

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

cunits = '[mm/step]';
uu = uu*1000; % m -> mm
vv = vv*1000;
mm = mm*1000;
ulim = ulim*1000;
vlim = vlim*1000;
mlim = mlim*1000;

% create figures and axes
[ax_top, ax_mid, ax_bot] = create_figure();

% plot displacement magnitude and direction
axes(ax_top)
imagesc(xx, yy, mm, 'AlphaData', ~isnan(mm)); 
format_axes(gca, xlim, xunits, 0, ylim, yunits, 1, mlim, cunits, ...
    sprintf('displacement vector magnitude, step = %.1f', step));

% plot x-direction displacement magnitude
axes(ax_mid);
imagesc(xx, yy, uu, 'AlphaData', ~isnan(uu)); 
format_axes(gca, xlim, xunits, 0, ylim, yunits, 1, ulim, cunits, ...
    sprintf('x-direction displacement component, step = %.1f', step));

% plot y-direction displacement magnitude
axes(ax_bot);
imagesc(xx, yy, vv, 'AlphaData', ~isnan(vv)); 
format_axes(gca, xlim, xunits, 1, ylim, yunits, 1, vlim, cunits, ...
    sprintf('y-direction displacement component, step = %.1f', step));

% debug: add quiver plot to top
axes(ax_top)


keyboard

% generate quiver grid: staggered regular grid within the roi
q_bnd_frac = 0.05;
nxq = 10; 
nyq = 5; 
xq = linspace(xlim(1)+range(xlim)*q_bnd_frac, xlim(2)-range(xlim)*q_bnd_frac, nxq);
yq = linspace(ylim(1)+range(ylim)*q_bnd_frac, ylim(2)-range(ylim)*q_bnd_frac, nyq);
[xq, yq] = meshgrid(xq, yq);
[xgrid, ygrid] = meshgrid(xx, yy);
roi = ~isnan(uu);
roiq = interp2(xgrid, ygrid, roi, xq, yq, 'nearest');
xq = xq(roiq);
yq = yq(roiq);

% interpolate vectors and plot
uq = spline2d(xq, yq, xgrid(roi), ygrid(roi), uu(roi), tension);
vq = spline2d(xq, yq, xgrid(roi), ygrid(roi), vv(roi), tension);

hold on
quiver(xq, yq, uq, vq);



% % generate quiver plot

% 
% 
% [qc, qr] = meshgrid(1:q_spc_x:nx, 1:q_spc_y:ny);
% qr(:,2:2:end) = qr(:,2:2:end)+q_spc_y/2; % stagger
% qk = sub2ind(size(uu), qr, qc);
% 
% roi = ~isnan(uu);
% qk = qk(roi(qk)); 
% 
% figure
% plot(qc(:), qr(:), '.k');

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

function [] = format_axes(ax, xlim, xunits, xshow, ylim, yunits, yshow, ...
                climits, cunits, title_str)

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

title(title_str);

caxis(climits);
h = colorbar;
ylabel(h, cunits);
    
end
