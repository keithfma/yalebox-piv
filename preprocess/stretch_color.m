function [img_stretch, stretch_min, stretch_rng] = stretch_color(img, mask, upper_limit, lower_limit)
%
% Compute and/or apply stretch correction to make colors use the full
% available range. Stretch correction vectors as a function of x are
% computed from a LOESS curve fit to per-column quantiles. The stretched
% image is (img-stretch_min)/stretch_rng, truncated to the range [0, 1].
% The stretch is computed independently for each color band.
%
%
% Usage:
%   Compute and apply:
%   Compute only:
%   Correct only:
%
% Arguments:
%   img = 
%   mask = 
%   upper_limit = 
%   lower_limit = 
%
% Keith Ma, July 2015

% DEBUG: grayscale
img = rgb2gray(img);

%% get correction from smoothed per-column quantiles

ncol = size(img,2);
x = 1:ncol;
lower = nan(size(x));
upper = nan(size(x));
for i = x
    img_sub = img(:,i);
    mask_sub = mask(:,i);
    tmp = quantile(img_sub(mask_sub), [lower_limit, upper_limit]);
    lower(i) = tmp(1);
    upper(i) = tmp(2);
end

x_fit = downsample(x,10);
if x_fit(end) ~= x(end); 
    x_fit(end+1) = x(end);
end

lower_fit = loess(x, lower, x_fit, 0.25, 1, 0);
upper_fit = loess(x, upper, x_fit, 0.25, 1, 0);

lower_smooth = interp1(x_fit, lower_fit, x, 'linear');
upper_smooth = interp1(x_fit, upper_fit, x, 'linear');

stretch_min = lower_smooth;
stretch_rng = upper_smooth-lower_smooth;

%% apply stretch

img_stretch = img-repmat(stretch_min, size(img,1), 1);
img_stretch = img_stretch./repmat(stretch_rng, size(img,1), 1);
img_stretch(:) = min(1, max(0, img_stretch(:)));


