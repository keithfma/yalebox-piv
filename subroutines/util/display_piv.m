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

if nargin < 5
    prefix = '';
end

% extract key variables from the PIV results, apply ROI mask
xx = piv_result.x;
yy = piv_result.y;
uu = piv_result.u;
vv = piv_result.v;
roi = ~isnan(uu);

% get data limits
xlim = [min(xx(roi)), max(xx(roi))] + range(xx(roi))*[-0.05, 0.05];
ylim = [min(yy(roi)), max(yy(roi))] + range(yy(roi))*[-0.15, 0.15];

% plot displacement magnitude and direction, [mm/step]
figure
ax1 = subplot(3,1,1);
mm = sqrt(uu.^2 + vv.^2);
imagesc(xx(1,:), yy(:,1), mm*1000, 'AlphaData', ~isnan(mm)); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
title([prefix, 'Displacement Magnitude']);

% plot x-direction displacement magnitude, [mm/step]
ax2 = subplot(3,1,2);
imagesc(xx(1,:), yy(:,1), uu*1000, 'AlphaData', ~isnan(uu)); 
colormap(gca, flipud(colormap)); % flow is in negative x direction
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
title([prefix, 'X-Direction Displacement']);

% plot y-direction displacement magnitude, [mm/step]
ax3 = subplot(3,1,3);
imagesc(xx(1,:), yy(:,1), vv*1000, 'AlphaData', ~isnan(vv)); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
title([prefix, 'Y-Direction Displacement']);

linkaxes([ax1, ax2, ax3]);
