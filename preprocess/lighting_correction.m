function [imgc] = lighting_correction(img, mask, width, step, quant, lfrac)
%
% Compute and apply color baseline and scale corrections for masked
% sandbox image.
%
% Arguments:
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

%% windowed quantiles


%% apply corrections

%% debug
imgc = [];
fprintf('width = %.1f\n', width);
fprintf('step = %.1f\n', step);
fprintf('quant = [%.3f, %.3f]\n', quant(1), quant(2));
fprintf('lfrac = %.3f\n', lfrac);
