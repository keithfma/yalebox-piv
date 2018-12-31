function display_image_pair(xx, yy, ini, ini_roi, fin, fin_roi)
% function display_image_pair(xx, yy, ini, ini_roi, fin, fin_roi)
%
% Simple plot of image pair with area outside ROI set to 0
%
% Arguments:
%   xx, yy: Vector, world coordinate vectors for image matrices
%   ini, fin: 2D matrix, intensity values for initial and final images 
%   ini_roi, fin_roi: 2D matrix, logical flags for inital and final images,
%       indicating if a pixel is sand (1) or not (0)
% %

figure 

ini_masked = ini;
ini_masked(~ini_roi) = 0;

ax1 = subplot(2,1,1);
imagesc(xx, yy, ini_masked);
colormap('gray');
ax1.YDir = 'normal';
axis equal tight
title('Initial Image');

fin_masked = fin;
fin_masked(~fin_roi) = 0;

ax2 = subplot(2,1,2);
imagesc(xx, yy, fin);
colormap('gray');
ax2.YDir = 'normal';
axis equal tight
title('Final Image');

linkaxes([ax1, ax2]);
