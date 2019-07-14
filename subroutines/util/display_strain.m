function display_strain(strain_result, piv_mask, prefix)
% function display_strain(strain_result, piv_mask, prefix)
%
% Simple plot of strain parameters derived from PIV displacement fields,
% intended for rapid inspection of results
%
% Arguments:
%   strain_result: Struct, results for strain analysis, as returned by PIV analysis as returned by piv()
%   prefix: String, optional prefix for subplot titles
% %

% FIXME: the piv_mask should be burned into the strain results as well

if nargin < 3
    prefix = '';
end

% apply mask
Dd = strain_result.Dd; Dd(~piv_mask) = NaN;
Dv = strain_result.Dv; Dv(~piv_mask) = NaN;
spin = strain_result.spin; spin(~piv_mask) = NaN;


% get data limits
[xx, yy] = meshgrid(strain_result.x, strain_result.y);
roi = ~isnan(strain_result.Dd);
xlim = [min(xx(roi)), max(xx(roi))];
ylim = [min(yy(roi)), max(yy(roi))];
Ddlim = prctile(strain_result.Dd(:), [5, 95]);
Dvlim = prctile(strain_result.Dv(:), [5, 95]);
spinlim = prctile(strain_result.spin(:), [5, 95]);

% plot Dd
figure

a1 = subplot(3,1,1);
imagesc(strain_result.x, strain_result.y, Dd, ...
    'AlphaData', ~isnan(Dd)); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
caxis(Ddlim);
title([prefix, 'Dd']);

% plot Dv
a2 = subplot(3,1,2);
imagesc(strain_result.x, strain_result.y, Dv, ...
    'AlphaData', ~isnan(Dv)); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
caxis(Dvlim);
title([prefix, 'Dv']);

% plot spin
a3 = subplot(3,1,3);
imagesc(strain_result.x, strain_result.y, spin, ...
    'AlphaData', ~isnan(spin)); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
caxis(spinlim);
title([prefix, 'spin']);

linkaxes([a1, a2, a3]);
