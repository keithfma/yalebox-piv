function [out] = correct_apply(rgb, mask, baseline, scale, tiny)
%
% Apply mask and color corrections computed using the correct_mask() and
% correct_color() functions.
%
% Arguments:
%
%   rgb = 3D matrix, uint8, a 24-bit "Truecolor" RGB image, as read into
%       MATLAB with imread()
%
%   mask = 2D matrix, double, TRUE where there is sand, FALSE elsewhere
%
%   baseline = 2D matrix, size(img, 1) x size(img, 2), baseline to be subtracted
%       from the image
%
%   tiny = Scalar, double, minimum intensity value for sand pixels, should
%       be greater than 0 to differentiate from mask pixels, default = 1e-3
%
%   scale = 2D matrix, size(img, 1) x size(img, 2), scaling factor to multiply 
%       the image by
%
%   out = 2D matrix, double, corrected image intensity, normalized to the range
%       [0, 1]
%
% Keith Ma, July 2015

% set defaults
if isempty(tiny); tiny = 1e-3; end

% check for sane inputs
assert(isa(rgb, 'uint8') && size(rgb,3) == 3, ...
     'rgb is not a Truecolor (24-bit) RGB image');
assert(isa(mask, 'logical'), ...
    'mask is not of type logical ');
assert(size(mask,1) == size(rgb, 1) && size(mask,2) == size(rgb, 2), ...
    'mask and rgb dimensions do not match');
assert(isa(baseline, 'double'), 'baseline is not of type double');
assert(size(baseline,1) == size(rgb, 1) && size(baseline,2) == size(rgb, 2), ...
    'baseline and rgb dimensions do not match');
assert(isa(scale, 'double'), 'scale is not of type double');
assert(size(scale,1) == size(rgb, 1) && size(scale,2) == size(rgb, 2), ...
    'scale and rgb dimensions do not match');
assert(isa(scale, 'double'), 'scale is not of type double');
assert(isa(tiny, 'double') && numel(tiny) == 1, ...
    'tiny is not a scalar double');

% apply corrections
hsv = rgb2hsv(rgb);
out = hsv(:,:,3);
out = out-baseline;
out = out.*scale;
out(out<tiny) = tiny;
out(~mask) = 0;