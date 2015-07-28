function [baseline, scale] = correct_color(rgb, mask, width, npts, quant, lfrac)
% function [baseline, scale] = correct_color(rgb, mask, width, npts, quant, lfrac)
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
%   rgb = 3D matrix, double, TRUECOLOR RGB 
%
%   mask = 2D matrix, double, TRUE where there is sand, FALSE elsewhere
%
%   width = Scalar, double, width of the sliding window for quantile
%       calculations, in (whole) pixels
%
%   npts = Scalar, integer, number of linearly spaced points at which
%       quantiles are computed.
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
if isempty(npts); npts = 100; end
if isempty(quant); quant = [0.01, 0.99]; end
if isempty(lfrac); lfrac = 0.2; end

% check for sane inputs
assert(isa(rgb, 'uint8') && size(rgb,3) == 3, ...
     'rgb is not a Truecolor (24-bit) RGB image');
assert(isa(mask, 'logical'), ...
    'mask is not of type logical ');
assert(size(mask,1) == size(rgb, 1) && size(mask,2) == size(rgb, 2), ...
    'mask and rgb dimensions do not match');
assert(numel(width) == 1 && mod(width,1) == 0, ...
    'width is not a scalar integer');
assert(numel(npts) == 1 && mod(npts,1) == 0, ...
    'npts is not a scalar integer');
assert(numel(quant) == 2 && max(quant) >=  0 && min(quant) <= 1, ...
    'quant is not a 2-element vector in the range 0-1');
assert(numel(lfrac) == 1 && max(lfrac) >=  0 && min(lfrac) <= 1, ...
    'lfrac is not a scalar in the range 0-1');

% convert image to 2D intensity matrix
hsv = rgb2hsv(rgb);
val = hsv(:,:,3);
[nrow, ncol] = size(val);

% create bins, truncated at image edges

x = linspace(1, ncol, npts);
x0 = max(1, round(x-width));
x1 = min(ncol, round(x+width));

% windowed quantiles for intensity
top = nan(1, npts);
bot = nan(1, npts);
for i = 1:npts
    val_win = val(:,x0(i):x1(i));
    mask_win = mask(:,x0(i):x1(i));
    tmp = quantile(val_win(mask_win), quant);
    bot(i) = tmp(1);
    top(i) = tmp(2);
end

% smooth and upscale quantiles
top = loess(x, top, x, lfrac, 1, 0);
top = interp1(x, top, 1:ncol, 'linear');

bot = loess(x, bot, x, lfrac, 1, 0);
bot = interp1(x, bot, 1:ncol, 'linear');

% correction matrices from smooth quantiles
baseline = repmat(bot,          nrow, 1);
scale =    repmat(1./(top-bot), nrow, 1);