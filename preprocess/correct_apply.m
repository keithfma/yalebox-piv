function [imgc] = correct_apply(img, mask, baseline, scale)
%
% Apply mask and color corrections computed using the correct_mask() and
% correct_color() functions.
%
% Arguments:
%
%   img = 3D matrix, double, RGB image, normalized to the range [0, 1]
%   mask = 2D matrix, double, TRUE where there is sand, FALSE elsewhere
%   baseline = Vector, 1x(size(img,2)), baseline to be subtracted from the image
%   scale = Vector, 1x(size(img,2)), scaling factor to multiply the image by
%   imgc = 3D matrix, double, corrected RGB image, normalized to the range
%       [0, 1]
%
% Keith Ma, July 2015

%% check for sane inputs

assert(isa(img, 'double'), 'img is not of type double');
assert(size(img,3) == 3, 'img is not a 3D array (RGB)');

assert(isa(mask, 'logical'), 'mask is not of type logical');
assert(size(mask,1) == size(img, 1) && size(mask,2) == size(img, 2), ...
    'mask and image dimensions do not match');

assert(isa(baseline, 'double'), 'baseline is not of type double');
assert(size(baseline, 2) == size(img, 2), ...
    'baseline does not match the number of columns in img');

assert(isa(scale, 'double'), 'scale is not of type double');
assert(size(scale, 2) == size(img, 2), ...
    'scale does not match the number of columns in img');


%% apply corrections

imgc = img;
imgc(repmat(~mask, [1, 1, 3])) = 0;
imgc = bsxfun(@minus, imgc, baseline);
imgc = bsxfun(@times, imgc, scale);
imgc = min(1, max(0, imgc));
