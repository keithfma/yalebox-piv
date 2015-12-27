function [] = piv_test_util_plot(xx, yy, uu, vv, displ, Dd)
%
% Simple plots of select PIV outputs
%
% %

figure
imagesc(xx, yy, displ);
colorbar
hold on
quiver(xx, yy, uu, vv, 2, '-k');
axis equal
axis tight
hold off
title('Displacement')

figure
imagesc(xx, yy, Dd)
colorbar
axis equal
axis tight
title('Strain (2nd Invariant)')