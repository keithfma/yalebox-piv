function [] = test_piv_util_plot(xx, yy, uu, vv, displ, Dd)
%
% Simple plots of select PIV outputs
%
% %
fig_position = [0.4312    0.4154    0.3558    0.4180];

figure('Units', 'Normalized', 'Position', fig_position);
imagesc(xx, yy, displ, 'AlphaData', ~isnan(displ));
set(gca, 'YDir', 'normal');
colorbar
hold on
quiver(xx, yy, uu, vv, 2, '-k');
axis equal
axis tight
hold off
title('Velocity Magnitude and Direction')

figure('Units', 'Normalized', 'Position', fig_position);
imagesc(xx, yy, Dd, 'AlphaData', ~isnan(Dd))
set(gca, 'YDir', 'normal');
colorbar
axis equal
axis tight
title('Strain (2nd Invariant)')