function plot_image_pair(xx, yy, ini, ini_roi, fin, fin_roi)
% function plot_image_pair(xx, yy, ini, ini_roi, fin, fin_roi)
%
% Simple plots of image pair, including one figure with images and a separate
% figure with the corresponding ROI masks
%
% Arguments:
%   xx, yy: Vector, world coordinate vectors for image matrices
%   ini, fin: 2D matrix, intensity values for initial and final images 
%   ini_roi, fin_roi: 2D matrix, logical flags for inital and final images,
%       indicating if a pixel is sand (1) or not (0)
% %


figure 

ax1 = subplot(2,1,1);
imagesc(xx, yy, ini);
colormap('gray');
ax1.YDir = 'normal';
axis equal tight
title('Initial Image');

ax2 = subplot(2,1,2);
imagesc(xx, yy, fin);
colormap('gray');
ax2.YDir = 'normal';
axis equal tight
title('Final Image');

linkaxes([ax1, ax2]);

figure 

ax3 = subplot(2,1,1);
imagesc(xx, yy, ini_roi);
colormap('gray');
ax3.YDir = 'normal';
axis equal tight
title('Initial Image ROI');

ax4 = subplot(2,1,2);
imagesc(xx, yy, fin_roi);
colormap('gray');
ax4.YDir = 'normal';
axis equal tight
title('Final Image ROI');

linkaxes([ax3, ax4]);