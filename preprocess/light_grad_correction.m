function [correct, img_correct] = light_grad_correction(img, mask, loess_frac, correct);
% function [correct, img_correct] = light_grad_correction(img, mask, loess_frac, correct);
%
% Compute and/or apply lighting correction factor as a function of x for undesirable
% x-direction lighting gradients. 
%
% Arguments:
%   img = 3D matrix, double. RGB image converted to double, range is arbitrary
%   mask = 1D matrix, logical. TRUE where pixels are sand, FALSE elsewhere,
%       as produced by sand_mask()
%   loess_frac = Scalar, double. Parameter to loess(), defined there as
%       "smoothing (fraction of series to include in weighted local fit)
%       typically 0.25 to 1.0"
%   correct = Vector, double. Additive correction factor to remove
%       gradients, can be both input and output.
%   img_correct = 3D matrix, double. Corrected RGB image.
%
% Usage:
%   Compute and apply:
%       [correct, img_correct] = light_grad_correction(img, mask, loess_frac);
%   Compute only:
%       [correct] = light_grad_correction(img, mask, loess_frac);
%   Apply only:
%       [~, img_correct] = light_grad_correction(img, mask, [], correct);
% 
% The correction is computed by (1) converting the image to normalized
% grayscale, (2) computing mean intensity of sand pixels as a function of
% x, (3) fitting a smooth line to the mean intensity data using loess, (4)
% computing the offset needed to flatten the smoothed intensity curve.
%
% The correction is applied as follows:
%   img_correct = img + repmat(correct, size(img,1), 1);
%   img_correct = min(img_max, max(img_min, img_corrected);
%
% Keith Ma, July 2015

%% parse mode

if nargin == 3 && ~isempty(loess_frac) && nargout == 2        
    get_correct = 1;
    apply_correct = 1;
    
elseif nargin == 3 && ~isempty(loess_frac) && nargout == 1
    get_correct = 1;
    apply_correct = 0;
    
elseif nargin == 4 && isempty(loess_frac) && nargout == 2
    get_correct = 0;
    apply_correct = 1;
    
else
    error('Incorrect usage. Type "help light_grad_correction" for more details');
end

%% check for sane inputs
    
assert(isa(img, 'double'), 'img is not of type double');
assert(size(img,3) == 3, 'img is not a 3D array (RGB)');

assert(isa(mask, 'logical'), 'mask is not of type logical');
assert(size(mask,1) == size(img, 1) && size(mask,2) == size(img, 2), ...
    'mask and image dimensions do not match');

if get_correct
    assert(numel(loess_frac) == 1, 'loess_frac is not a scalar');
    assert(isa(loess_frac, 'double'), 'loess_frac in not of type double');
end

if nargin == 4
    assert(min(size(correct)) == 1 && max(size(correct)) == size(img,2), ...
        'correct is not a vector with size(img,1) elements');
end

%% extract average intensity of sand pixels as (x, intensity) points

img_gray = rgb2gray(img);
img_gray(~mask) = NaN;
intens = nanmean(img_gray,1);
 
%% fit smooth model
% Note: higher-order interpolation differs by O(10^-6) from linear

x = 1:size(img,2);
x_fit = downsample(x,10);
if x_fit(end) ~= x(end); 
    x_fit(end+1) = x(end);
end
 
intens_fit = loess(x, intens, x_fit, loess_frac, 1, 0);
intens_smooth = interp1(x_fit, intens_fit, x, 'linear');

%% compute and/or apply correction

if get_correct
    mean_intens_smooth = mean(intens_smooth);
    correct = mean_intens_smooth - intens_smooth;    
end

if apply_correct
    img_correct = img + repmat(correct, size(img_gray,1), 1, 3);    
end
