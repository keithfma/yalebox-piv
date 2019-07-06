function display_piv(piv_result, prefix)
% function display_piv(piv_result, prefix)
%
% Simple plot of PIV displacement fields, intended for rapid inspection of
% results
%
% Arguments:
%   piv_result: Struct, results for PIV analysis as returned by piv()
%   prefix: String, optional prefix for subplot titles
% %

update_path('piv');

if nargin < 2
    prefix = '';
end

% extract key variables from the PIV results
xx = piv_result.x;
yy = piv_result.y;
uu = piv_result.u;
vv = piv_result.v;
qual = piv_result.quality;
mask = piv_result.mask;

% get data limits
xlim = [min(xx(:)), max(xx(:))];
ylim = [min(yy(:)), max(yy(:))];

hf = figure;

% plot displacement magnitude and direction, [mm/step]
% FIXME: ain't no direction here yet...
% FIXME: draw mask boundaries on these plots
ax1 = subplot(4, 1, 1);
mm = sqrt(uu.^2 + vv.^2);
imagesc(xx(1,:), yy(:,1), mm*1000); 
colorbar;
draw_mask(xx, yy, mask);
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
title([prefix, 'Displacement Magnitude']);

% plot x-direction displacement magnitude, [mm/step]
ax2 = subplot(4, 1, 2);
imagesc(xx(1,:), yy(:,1), uu*1000); 
colormap(gca, flipud(colormap)); % flow is in negative x direction
colorbar;
draw_mask(xx, yy, mask);
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
title([prefix, 'X-Direction Displacement']);

% plot y-direction displacement magnitude, [mm/step]
ax3 = subplot(4, 1, 3);
imagesc(xx(1,:), yy(:,1), vv*1000); 
colorbar;
draw_mask(xx, yy, mask);
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
title([prefix, 'Y-Direction Displacement']);

% plot quality flags, [categorical]
ax4 = subplot(4, 1, 4);
imagesc(xx(1,:), yy(:,1), Quality.to_uint8(qual), ...
    'AlphaData', qual ~= Quality.Valid);
[members, names] = enumeration('Quality');
names(members == Quality.Valid) = [];
members(members == Quality.Valid) = [];
values = Quality.to_uint8(members);
caxis([min(values), max(values)])
colormap(gca, lines(numel(members)));
hcb = colorbar;
hcb.Ticks = values;
hcb.TickLabels = names;
draw_mask(xx, yy, mask);
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
title([prefix, 'Quality Flags']);

% some final formatting
linkaxes([ax1, ax2, ax3, ax4]);
hf.Units = 'Normalized';
hf.OuterPosition = [0, 0, 1, 1];


function draw_mask(xx, yy, mask)

hold on
mask_bounds = bwboundaries(mask);
for ii = 1:numel(mask_bounds)
    mask_x = xx(1, mask_bounds{ii}(:, 2));
    mask_y = yy(mask_bounds{ii}(:, 1), 1);
    plot(mask_x, mask_y, '-k');
end