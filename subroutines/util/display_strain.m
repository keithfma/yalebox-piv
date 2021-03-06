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

% unpack data
Dd = strain_result.Dd;
Dv = strain_result.Dv;
spin = strain_result.spin;


% get data limits
[xx, yy] = meshgrid(strain_result.x, strain_result.y);
xlim = [min(xx(piv_mask)), max(xx(piv_mask))];
ylim = [min(yy(piv_mask)), max(yy(piv_mask))];
Ddlim = prctile(strain_result.Dd(piv_mask), [15, 85]);
Dvlim = prctile(strain_result.Dv(piv_mask), [15, 85]);
spinlim = prctile(strain_result.spin(piv_mask), [15, 85]);

hf = figure;

% plot Dd
a1 = subplot(3,1,1);
imagesc(strain_result.x, strain_result.y, Dd); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
caxis(Ddlim);
title([prefix, ' Dd'], 'Interpreter', 'none');

% plot Dv
a2 = subplot(3,1,2);
imagesc(strain_result.x, strain_result.y, Dv); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
caxis(Dvlim);
title([prefix, ' Dv'], 'Interpreter', 'none');

% plot spin
a3 = subplot(3,1,3);
imagesc(strain_result.x, strain_result.y, spin); 
colorbar;
axis equal tight
set(gca, 'YDir', 'normal', 'XLim', xlim, 'YLim', ylim);
caxis(spinlim);
title([prefix, ' spin'], 'Interpreter', 'none');

% some final formatting
linkaxes([a1, a2, a3]);
hf.Units = 'Normalized';
hf.OuterPosition = [0, 0, 1, 1];