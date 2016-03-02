function [] = plot_displacement(piv_file, index, xlimits, ylimits, ulimits, vlimits, mlimits)
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
% xlimits, ylimits, ulimits, vlimits, mlimits = Vector, length==2, [minimum, maximum] values for the x-axis, y-axis, u-displacement, v-displacement, and
%   displacement magnitude. Values beyond the range of the data will be
%   truncated to the data limits. Thus, to span the data, one could use [-inf, inf].
% %

% debug: hard-code input arguments
piv_file = '~/Documents/dissertation/yalebox-exp-fault/data/fault_ss_01/piv/fault_ss_01_sidef.displ.nc';
index = 100;
xlimits = [-inf, inf];
ylimits = [-inf, 0.1];
ulimits = [-inf, inf];
vlimits = [0, inf];
mlimits = [-inf, inf];

% check for sane inputs
validateattributes(piv_file, {'char'}, {'vector'});
validateattributes(index, {'numeric'}, {'scalar', 'positive', 'integer'});
validateattributes(xlimits, {'numeric'}, {'vector', 'numel', 2, 'increasing'});
validateattributes(ylimits, {'numeric'}, {'vector', 'numel', 2, 'increasing'});
validateattributes(ulimits, {'numeric'}, {'vector', 'numel', 2, 'increasing'});
validateattributes(vlimits, {'numeric'}, {'vector', 'numel', 2, 'increasing'});
validateattributes(mlimits, {'numeric'}, {'vector', 'numel', 2, 'increasing'});

% read data from netCDF
xx = ncread(piv_file, 'x'); nx = numel(xx);
yy = ncread(piv_file, 'y'); ny = numel(yy);
step = ncread(piv_file, 'step', index, 1); % unused
uu = squeeze(ncread(piv_file, 'u', [1, 1, index], [nx, ny, 1]))';
vv = squeeze(ncread(piv_file, 'v', [1, 1, index], [nx, ny, 1]))';
mm = sqrt(uu.*uu+vv.*vv);

% truncate axis limits
xlimits(1) = max(xlimits(1), min(xx)); xlimits(2) = min(xlimits(2), max(xx));
ylimits(1) = max(ylimits(1), min(yy)); ylimits(2) = min(ylimits(2), max(yy));
ulimits(1) = max(ulimits(1), min(uu(:))); ulimits(2) = min(ulimits(2), max(uu(:)));
vlimits(1) = max(vlimits(1), min(vv(:))); vlimits(2) = min(vlimits(2), max(vv(:)));
mlimits(1) = max(mlimits(1), min(mm(:))); mlimits(2) = min(mlimits(2), max(mm(:)));

% convert units
xunits = '[cm]';
xx = xx*100; % m -> cm
xlimits = xlimits*100;

yunits = '[cm]';
yy = yy*100; % m -> cm
ylimits = ylimits*100;

cunits = '[mm/step]';
uu = uu*1000; % m -> mm
vv = vv*1000;
mm = mm*1000;
ulimits = ulimits*1000;
vlimits = vlimits*1000;
mlimits = mlimits*1000;

% create figures and axes
[ax_top, ax_mid, ax_bot] = create_figure();

% plot displacement magnitude and direction
axes(ax_top)
imagesc(xx, yy, mm, 'AlphaData', ~isnan(mm)); 
format_axes(gca, xlimits, xunits, 0, ylimits, yunits, 1, mlimits, cunits, ...
    sprintf('displacement vector magnitude, step = %.1f', step));

% plot x-direction displacement magnitude
axes(ax_mid);
imagesc(xx, yy, uu, 'AlphaData', ~isnan(uu)); 
format_axes(gca, xlimits, xunits, 0, ylimits, yunits, 1, ulimits, cunits, ...
    sprintf('x-direction displacement component, step = %.1f', step));

% plot y-direction displacement magnitude
axes(ax_bot);
imagesc(xx, yy, vv, 'AlphaData', ~isnan(vv)); 
format_axes(gca, xlimits, xunits, 1, ylimits, yunits, 1, vlimits, cunits, ...
    sprintf('y-direction displacement component, step = %.1f', step));

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

function [] = format_axes(ax, xlimits, xunits, xshow, ylimits, yunits, yshow, ...
                climits, cunits, title_str)

axis equal
grid on

set(ax, 'XLim', xlimits, ...
        'YLim', ylimits, ...
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
