function [] = test_piv_util_plot_error(uerr, verr)
%
% Display a quick plot of the error matrices
 
figure('units', 'normalized', 'position', [0.05, 0.3, 0.9, 0.5]);

subplot(1, 3, 1)
imagesc(uerr);
set(gca, 'YDir', 'normal');
colorbar
title('uu error')

subplot(1, 3, 2)
imagesc(verr);
set(gca, 'YDir', 'normal');
colorbar
title('vv error')

subplot(1, 3, 3)
imagesc(sqrt(uerr.^2+verr.^2));
set(gca, 'YDir', 'normal');
colorbar
title('error magnitude')

end