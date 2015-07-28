function [baseline, scale] = correct_color(img, mask, width, step, quant, lfrac)
% function [baseline, scale] = correct_color(img, mask, width, step, quant, lfrac)
%
% Compute color baseline and scale corrections for masked sandbox image.
% Specifically, sandbox images frequently have lateral gradients (with
% respect to x) in brightness and contrast. This function computes upper
% and lower quantiles as a function of x using a sliding window, and uses
% smooth curves fit to these data to adjust the baseline and range of the
% intensity data. The mask matrix makes it possible to include only the
% sand region in the quantiles. In pseudocode, the color correction is:
%
%   corrected_image = (original_image-baseline)*scale
%   corrected_image = min(1, corrected_image)
%   corrected_image = max(0, corrected_image)
%
% All arguments are required - default values are used where inputs are []
%
% Arguments:
%
%   img = 3D matrix, double, RGB image, normalized to the range [0, 1]
%
%   mask = 2D matrix, double, TRUE where there is sand, FALSE elsewhere
%
%   width = Scalar, double, width of the sliding window for quantile
%       calculations, in (whole) pixels
%
%   step = Scalar, double, spacing of windows for quantile calculations, in
%       (whole) pixels
%
%   quant = 2-element vector, double, [lower, upper] quantiles to be
%       calculated (and truncated by the correction), in the range 0 to 1
%
%   lfrac = Scalar, double, parameter to LOESS smooth curve fitting
%       algorithm, the fraction of the dataset to include in each local fit
%
%   baseline = Vector, 1x(size(img,2)), baseline to be subtracted from the image
%
%   scale = Vector, 1x(size(img,2)), scaling factor to multiply the image by
%
% Keith Ma, July 2015

% set defaults
if isempty(width); width = 200; end 
if isempty(step); step = 50; end
if isempty(quant); quant = [0.01, 0.99]; end
if isempty(lfrac); lfrac = 0.2; end

% check for sane inputs
assert(isa(img, 'double'), ...
    'img is not of type double');
assert(size(img,3) == 3, ...
    'img is not a 3D array (RGB)');
assert(isa(mask, 'logical'), ...
    'mask is not of type logical');
assert(size(mask,1) == size(img, 1) && size(mask,2) == size(img, 2), ...
    'mask and image dimensions do not match');
assert(numel(width) == 1, ...
    'width is not a scalar');
assert(mod(width,1) == 0, ...
    'width is not an integer');
assert(numel(step) == 1, ...
    'step is not a scalar');
assert(mod(step,1) == 0, ...
    'step is not an integer');
assert(numel(quant) == 2, ...
    'quant is not a 2-element vector');
assert(max(quant) >=  0 && min(quant) <= 1, ...
    'quant is not in the range 0-1');
assert(numel(lfrac) == 1, ...
    'lfrac is not a scalar');
assert(max(lfrac) >=  0 && min(lfrac) <= 1, ...
    'lfrac is not in the range 0-1');

%% get corrections from windowed quantiles

% convert image to grayscale
imgg = rgb2gray(img);

% get windows
ncol = size(img,2);
x = 1:step:ncol;
x0 = max(1, x-ceil(width/2));
x1 = min(ncol, x+ceil(width/2));

% compute quantiles
nx = numel(x);
top = nan(1, nx);
bot = nan(1, nx);
for i = 1:nx
    imgg_win = imgg(:,x0(i):x1(i));
    mask_win = mask(:,x0(i):x1(i));
    tmp = quantile(imgg_win(mask_win), quant);
    bot(i) = tmp(1);
    top(i) = tmp(2);
end

% smooth quantiles
top = loess(x, top, x, lfrac, 1, 0);
bot = loess(x, bot, x, lfrac, 1, 0);

% upscale smooth quantiles
top = interp1(x, top, 1:ncol, 'linear');
bot = interp1(x, bot, 1:ncol, 'linear');

% corrections
baseline = bot;
scale = 1./(top-bot);
