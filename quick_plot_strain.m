function quick_plot_strain(strain_result, prefix)
% function quick_plot_strain(strain_result, prefix)
%
% Simple plot of strain parameters derived from PIV displacement fields,
% intended for rapid inspection of results
%
% Arguments:
%   strain_result: Struct, results for strain analysis, as returned by PIV analysis as returned by piv()
%   prefix: String, optional prefix for subplot titles
% %

if nargin < 5
    prefix = '';
end

% get data limits
[xx, yy] = meshgrid(strain_result.x, strain_result.y);
roi = ~isnan(strain_result.Dd);
xlim = [min(xx(roi)), max(xx(roi))] + range(xx(roi))*[-0.05, 0.05];
ylim = [min(yy(roi)), max(yy(roi))] + range(yy(roi))*[-0.15, 0.15];

% plot Dd
figure

a1 = subplot(3,1,1);
imagesc(strain_result.x, strain_result.y, strain_result.Dd, ...
    'AlphaData', ~isnan(strain_result.Dd)); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
title([prefix, 'Dd']);

% plot Dv
a2 = subplot(3,1,2);
imagesc(strain_result.x, strain_result.y, strain_result.Dv, ...
    'AlphaData', ~isnan(strain_result.Dv)); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
title([prefix, 'Dv']);

% plot spin
a3 = subplot(3,1,3);
imagesc(strain_result.x, strain_result.y, strain_result.spin, ...
    'AlphaData', ~isnan(strain_result.spin)); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
title([prefix, 'spin']);

linkaxes([a1, a2, a3]);