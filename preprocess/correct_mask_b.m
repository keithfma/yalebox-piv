%function [mask] = correct_mask(image_rgb, hue_range, value_range, entropy_range,
%                               median_window, entropy_window, opening_radius)
%
% Create a logical mask for a color image that is TRUE where there is sand and
% FALSE elsewhere. This can be used to remove (set to 0) the background in a
% image prior to PIV analysis or other applications.
%
% Sand is identified by remapping to HSV colorspace, filtering and thresholding
% the "hue" and "value" bands, . Default values for the various parameters were
% effective for the test data used in development, but may well be inadequate
% for other images.
%
% Arguments:
%
%   im_rgb = 3D matrix, uint8, a 24-bit "Truecolor" RGB image, as read into
%       MATLAB with imread()
%
%   hue_range = 2-element vector, double, range [0, 1]. [minimum, maximum] HSV
%     "hue" included as sand in the mask. Default = [0.0, 0.5]
%
%   value_range = 2-element vector, double, range [0,1]. [minimum, maximum] HSV
%     "value" included as sand in the mask, . Default = [0.0, 0.5]
%
%   entropy_range = 2-element vector, double, range [0, 1]. [minimum, maximum]
%     entropy included as sand in the mask. Default = [0.5, 1.0]
%
%   median_window = scalar, integer, window size in pixels for median filter
%     applied to all bands
%
%   entropy_window = scalar, integer, window size in pixels for entropy filter
%
%   opening_radius = scalar, double, radius of disk structuring element used in
%     mophological opening filter
%
% Keith Ma, July 2015

% DEBUG
median_size = 25;
entropy_size = 9;
hue_range = [0, 0.4]; % min, max
value_range = [0, 0.5];
entropy_range = [0.5, 1]; 
morph_radius = 10;
% END DEBUG

% check for sane arguments, set default values

assert(isa(img_rgb, 'uint8'), 'img_rgb is not type "uint8"');
assert(size(img_rgb,3) == 3, 'img_rgb is not 3-dimensional');

% get hue, value, and entropy, normalized to the range [0, 1]
img_hsv = rgb2hsv(img_rgb); 
hue = img_hsv(:,:,1);
value = img_hsv(:,:,3);
entropy = entropyfilt(value, true(entropy_size));

% smooth, preseving edges with median filter
hue_m = medfilt2(hue, median_size*[1, 1]);
value_m = medfilt2(value, median_size*[1, 1]);
entropy_m = medfilt2(entropy, median_size*[1, 1]);

% normalize range to [0,1]
hue_m = hue_m-min(hue_m(:)); hue_m = hue_m./max(hue_m(:));
value_m = value_m-min(value_m(:)); value_m = value_m./max(value_m(:));
entropy_m = entropy_m-min(entropy_m(:)); entropy_m = entropy_m./max(entropy_m(:));

% threshold bands
hue_t = hue_m >= hue_range(1) & hue_m <= hue_range(2);
value_t = value_m >= value_range(1) & value_m <= value_range(2);
entropy_t = entropy_m >= entropy_range(1) & entropy_m <= entropy_range(2);

% create mask
mask = hue_t & value_t & entropy_t;

% fill holes - top and bottom
wall = true(1, size(mask,2));
mask_c = [wall; mask; wall];
mask_c = imfill(mask_c, 'holes');
mask_c = mask_c(2:end-1, :);

% fill holes - left and right
wall = true(size(mask, 1), 1);
mask_c = [wall, mask_c, wall];
mask_c = imfill(mask_c, 'holes');
mask_c = mask_c(:, 2:end-1);

% clean up edges with morphological opening filter 
disk = strel('disk', morph_radius);
mask_c = imopen(mask_c, disk);
