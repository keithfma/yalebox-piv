function [] = test_piv_util_plot_error(uerr, verr)
%
% Display a quick plot of the error matrices

fig_position_map = [0.4312    0.4154    0.3558    0.4180];
fig_position_hist = [0.4268    0.3672    0.3133    0.4479];
num_bins = 50;

figure('units', 'normalized', 'position', fig_position_map);
imagesc(uerr, 'AlphaData', ~isnan(uerr));
ax = gca;
ax.YDir= 'normal';
ax.XTickLabels = [];
ax.YTickLabels = [];
axis equal
axis tight
colorbar;
title('U Error [pixels]')

figure('units', 'normalized', 'position', fig_position_hist);
hist(uerr(~isnan(uerr)), num_bins);
ax = gca;
hh = findobj(gca,'Type','patch');
hh.LineStyle = 'none';
ax.XLim = [min(uerr(:)), max(uerr(:))];
title('U Error')
xlabel('Error [pixels]');
ylabel('Count [pixels]');

figure('units', 'normalized', 'position', fig_position_map);
imagesc(verr, 'AlphaData', ~isnan(verr));
ax = gca;
ax.YDir= 'normal';
ax.XTickLabels = [];
ax.YTickLabels = [];
axis equal
axis tight
colorbar
title('V Error [pixels]')

figure('units', 'normalized', 'position', fig_position_hist);
hist(verr(~isnan(verr)), num_bins);
ax = gca;
hh = findobj(gca,'Type','patch');
hh.LineStyle = 'none';
ax.XLim = [min(verr(:)), max(verr(:))];
title('V Error')
xlabel('Error [pixels]');
ylabel('Count [pixels]');

merr = sqrt(uerr.^2+verr.^2);

figure('units', 'normalized', 'position', fig_position_map);
imagesc(merr, 'AlphaData', ~isnan(verr) & ~isnan(merr));
ax = gca;
ax.YDir= 'normal';
ax.XTickLabels = [];
ax.YTickLabels = [];
axis equal
axis tight
colorbar
title('Velocity Magnitude Error [pixels]')

figure('units', 'normalized', 'position', fig_position_hist);
hist(merr(~isnan(merr)), num_bins);
ax = gca;
hh = findobj(gca,'Type','patch');
hh.LineStyle = 'none';
ax.XLim = [min(merr(:)), max(merr(:))];
title('Velocity Magnitude Error')
xlabel('Error [pixels]');
ylabel('Count [pixels]');

end