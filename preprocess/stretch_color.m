function [] = stretch_color(img, mask)
% Stretch colors to make full use of the available range.

upper_limit = 0.95;
lower_limit = 0.05;
box_width = 100; % pixels

img = rgb2gray(img);



% left-hand side
img_sub = img(:, 1:box_width/2);
mask_sub = mask(:, 1:box_width/2);
limits = quantile(img_sub(mask_sub), [lower_limit, upper_limit]);




