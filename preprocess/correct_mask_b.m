%function [mask] = correct_mask(rgb, hue_lim, value_lim, entropy_lim,
%                               median_window, entropy_window, opening_radius)
%
% Create a logical mask for a color image that is TRUE where there is sand and
% FALSE elsewhere. This can be used to remove (set to 0) the background in a
% image prior to PIV analysis or other applications. Sand is identified by
% remapping to HSV colorspace, filtering and thresholding the "hue" and "value"
% bands.
%
% All arguments MUST be provided. Default values will be used for any arguments
% set to the empty matrix []. The default values were effective for test data
% used in development, but are unlikely to work for all images. 
%
% Arguments:
%
%   im = 3D matrix, uint8, a 24-bit "Truecolor" RGB image, as read into
%       MATLAB with imread()
%
%   hue_lim = 2-element vector, double, range [0, 1]. [minimum, maximum] HSV
%     "hue" included as sand in the mask. Default = [0.0, 0.5]
%
%   value_lim = 2-element vector, double, range [0,1]. [minimum, maximum] HSV
%     "value" included as sand in the mask, . Default = [0.0, 0.5]
%
%   entropy_lim = 2-element vector, double, range [0, 1]. [minimum, maximum]
%     entropy included as sand in the mask. Default = [0.5, 1.0]
%
%   median_window = scalar, integer, window size in pixels for median filter
%     applied to all bands. Default = 25.
%
%   entropy_window = scalar, integer, window size in pixels for entropy filter.
%     Default = 9.
%
%   opening_radius = scalar, double, radius of disk structuring element used in
%     mophological opening filter. Default = 10.
%
% Keith Ma, July 2015

%% DEBUG
%median_size = 25;
%entropy_size = 9;
%hue_range = [0, 0.4]; % min, max
%value_range = [0, 0.5];
%entropy_range = [0.5, 1]; 
%morph_radius = 10;
%% END DEBUG

% set default values
if isempty(hue_lim); hue_lim = [0.0, 0.5]; end
if isempty(value_lim); value_lim = [0.0, 0.5]; end
if isempty(entropy_lim); entropy_lim = [0.5, 1.0]; end
if isempty(median_window); median_window = 25; end
if isempty(entropy_window); entropy_window = 9; end
if isempty(opening_radius); opening_radius = 10; end

% check for sane arguments, set default values
assert(isa(rgb, 'uint8'), 'img_rgb is not type "uint8"');
assert(size(rgb,3) == 3, 'img_rgb is not 3-dimensional');
assert(numel(hue_lim)==2 && min(size(hue_lim)) == 1, ...
    'hue_lim is not a 2-element vector');
assert(min(hue_lim) >= 0 && max(hue_lim) <= 1, ...
    'hue_lim is not in the range [0,1]');
assert(numel(value_lim)==2 && min(size(value_lim)) == 1, ...
    'value_lim is not a 2-element vector');
assert(min(value_lim) >= 0 && max(value_lim) <= 1, ...
    'value_lim is not in the range [0,1]');
assert(numel(entropy_lim)==2 && min(size(entropy_lim)) == 1, ...
    'entropy_lim is not a 2-element vector');
assert(min(entropy_lim) >= 0 && max(entropy_lim) <= 1, ...
    'entropy_lim is not in the range [0,1]');
assert(numel(median_window) == 1 && median_window > 0 && ...
    round(median_window) == median_window, ...
    'median_window is not a positive integer');
assert(numel(entropy_window) == 1 && entropy_window > 0 && ...
    round(entropy_window) == entropy_window, ...
    'entropy_window is not a positive integer');
assert(numel(opening_radius) == 1 && opening_radius > 0, ...
    'opening_radius is not a positive scalar');

% get hue, value, and entropy, normalized to the range [0, 1]
hsv = rgb2hsv(rgb); 
hue = hsv(:,:,1);
value = hsv(:,:,3);
entropy = entropyfilt(value, true(entropy_window));

% smooth, preseving edges with median filter
hue = medfilt2(hue, median_window*[1, 1]);
value = medfilt2(value, median_window*[1, 1]);
entropy = medfilt2(entropy, median_window*[1, 1]);

% normalize range to [0,1]
hue = hue-min(hue(:)); hue = hue./max(hue(:));
value = value-min(value(:)); value = value./max(value(:));
entropy = entropy-min(entropy(:)); entropy = entropy./max(entropy(:));

% threshold bands
hue_mask = hue >= hue_lim(1) & hue <= hue_lim(2);
value_mask = value >= value_lim(1) & value <= value_lim(2);
entropy_mask = entropy >= entropy_lim(1) & entropy <= entropy_lim(2);

% create mask
mask = hue_mask & value_mask & entropy_mask;

% fill holes - top and bottom
wall = true(1, size(mask,2));
mask = [wall; mask; wall];
mask = imfill(mask, 'holes');
mask = mask(2:end-1, :);

% fill holes - left and right
wall = true(size(mask, 1), 1);
mask = [wall, mask, wall];
mask = imfill(mask, 'holes');
mask = mask(:, 2:end-1);

% clean up edges with morphological opening filter 
disk = strel('disk', morph_radius);
mask = imopen(mask, disk);
