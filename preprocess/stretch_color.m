function [img_stretch, stretch_min, stretch_rng] = stretch_color(img, mask, trim, nbin, loess_frac)
% [img_stretch, stretch_min, stretch_rng] = ...
%       stretch_color(img, mask, trim, nbin, loess_frac)
%
% Compute and/or apply stretch correction to make colors use the full
% available range. Stretch correction vectors as a function of x are
% computed from a LOESS curve fit to windowed quantiles. The stretched
% image is (img-stretch_min)/stretch_rng, truncated to the range [0, 1].
% The stretch is computed independently for each color band.
%
% Usage:
%   Compute and apply:
%   Compute only:
%   Correct only:
%
% Arguments:
%   img = 
%   mask = 
%   trim =
%   nbin = 
%   loess_frac = 
%   
% Keith Ma, July 2015


%% check for sane arguments

get_stretch = 1;
apply_stretch = 1;
img_gray = rgb2gray(img);

%% get correction from windowed quantiles

if get_stretch
    
    edge = round(linspace(1, size(img,2), nbin));
    mid = edge(1:end-1)+diff(edge)/2;
    lower = nan(1, nbin-1);
    upper = nan(1, nbin-1);
    
    for i = 1:nbin-1
        img_sub = img_gray(:,edge(i):edge(i+1));
        mask_sub = mask(:,edge(i):edge(i+1));
        tmp = quantile(img_sub(mask_sub), trim);
        lower(i) = tmp(1);
        upper(i) = tmp(2);
    end
    
    lower_fit = loess(mid, lower, mid, loess_frac, 1, 0);
    upper_fit = loess(mid, upper, mid, loess_frac, 1, 0);
    
    lower_smooth = interp1(mid, lower_fit, 1:size(img,2), 'linear');
    upper_smooth = interp1(mid, upper_fit, 1:size(img,2), 'linear');
    
    stretch_min = lower_smooth;
    stretch_rng = upper_smooth-lower_smooth;
    
end

%% apply stretch

if apply_stretch
    
    img_stretch = img-repmat(stretch_min, size(img,1), 1, 3);
    img_stretch = img_stretch./repmat(stretch_rng, size(img,1), 1, 3);
    img_stretch(:) = min(1, max(0, img_stretch(:)));
    
end
