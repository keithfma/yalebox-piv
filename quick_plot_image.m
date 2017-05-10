function quick_plot_image(xx, yy, img, roi, prefix)
% function quick_plot_image(xx, yy, img, roi, prefix)
%
% Simple plot of image and its ROI mask, intended for rapid inspection of
% results
%
% Arguments:
%   xx, yy: Vector, world coordinate vectors for image matrices
%   img: 2D matrix, image intensity values
%   roi: 2D matrix, logical flags indicating if a pixel is sand (1) or not (0)
%   prefix: String, optional prefix for subplot titles
% %

if nargin < 5
    prefix = '';
end

figure 

ax1 = subplot(2,1,1);
imagesc(xx, yy, img);
colormap('gray');
ax1.YDir = 'normal';
axis equal tight
title([prefix, 'Image']);

ax2 = subplot(2,1,2);
imagesc(xx, yy, roi);
colormap('gray');
ax2.YDir = 'normal';
axis equal tight
title([prefix, 'Image ROI']);

linkaxes([ax1, ax2]);