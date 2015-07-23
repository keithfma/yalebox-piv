function [baseline, scale] = correct_color(img, mask, width, step, quant, lfrac)
% function [baseline, scale] = correct_color(img, mask, width, step, quant, lfrac)
%
% Compute color baseline and scale corrections for masked sandbox image.
%
% Arguments:
%
% img = 
% mask = 
% width =
% step = 
% quant = 
% lfrac = 
% baseline = 
% scale = 
%
% Keith Ma, July 2015

%% check for sane inputs

assert(isa(img, 'double'), 'img is not of type double');
assert(size(img,3) == 3, 'img is not a 3D array (RGB)');

assert(isa(mask, 'logical'), 'mask is not of type logical');
assert(size(mask,1) == size(img, 1) && size(mask,2) == size(img, 2), ...
    'mask and image dimensions do not match');

if nargin <= 2 || isempty(width)
    width = 200; % default
else
    assert(numel(width) == 1, 'width is not a scalar');
    assert(mod(width,1) == 0, 'width is not an integer');
end

if nargin <=3 || isempty(step)
    step = 50; %default
else
    assert(numel(step) == 1, 'step is not a scalar');
    assert(mod(step,1) == 0, 'step is not an integer');
end

if nargin <=4 || isempty(quant)
    quant = [0.01, 0.99]; %default
else
    assert(numel(quant) == 2, 'quant is not a 2-element vector');
    assert(max(quant) >=  0 && min(quant) <= 1, 'quant is not in the range 0-1');
end

if nargin <= 5 || isempty(lfrac)
    lfrac = 0.2;
else
    assert(numel(lfrac) == 1, 'lfrac is not a scalar');
    assert(max(lfrac) >=  0 && min(lfrac) <= 1, 'lfrac is not in the range 0-1');
end    

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
