function mask = yalebox_prep_mask_auto(rgb, hue_lim, value_lim, entropy_lim, median_window, entropy_window, morph_radius, show)
%
% function mask = yalebox_prep_mask_auto(rgb, ...
%                         hue_lim, value_lim, entropy_lim, ...
%                         median_window, entropy_window, morph_radius, show)
%
% Create a logical mask for a color image that is TRUE where there is sand and
% FALSE elsewhere. This can be used to remove (set to 0) the background in a
% image prior to PIV analysis or other applications. Sand is identified by
% remapping to HSV colorspace, filtering and thresholding the "hue" and "value"
% bands.
%
% Arguments:
%
%   rgb = 3D matrix, uint8, a 24-bit "Truecolor" RGB image, as read into
%       MATLAB with imread()
%
%   hue_lim = 2-element vector, double, range [0, 1]. [minimum, maximum] HSV
%     "hue" included as sand in the mask.
%
%   value_lim = 2-element vector, double, range [0,1]. [minimum, maximum] HSV
%     "value" included as sand in the mask.
%
%   entropy_lim = 2-element vector, double, range [0, 1]. [minimum, maximum]
%     entropy included as sand in the mask. 
%
%   median_window = scalar, integer, window size in pixels for median filter
%     applied to mask.
%
%   entropy_window = scalar, integer, window size in pixels for entropy filter.
%
%   morph_radius = scalar, double, radius of disk structuring element used in
%     mophological opening/closing filter. 
%
%   show = Scalar, logical, set to 1 (true) to plot the mask bands, used to
%       facilitate the parameter selection process, default = false.
%
%   mask = 2D matrix, logical, true where there is sand and false
%       elsewhere.
%
% Keith Ma, July 2015

% check for sane arguments, set default values
%assert(isa(rgb, 'uint8') && size(rgb,3) == 3, ...
%    'img_rgb is not a 24-bit RGB image');
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
assert(numel(morph_radius) == 1 && morph_radius > 0, ...
    'morph_radius is not a positive scalar');
if nargin == 7; show = false; end
assert(numel(show) == 1 && ismember(show, [0, 1]), ...
    'show is not 1 or 0');

% get hue, value, and entropy, normalized to the range [0, 1]
hsv = rgb2hsv(rgb); 
hue = hsv(:,:,1);
value = hsv(:,:,3);
entropy = entropyfilt(value, true(entropy_window));

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

% smooth, preseving edges with median filter
mask = medfilt2(mask, median_window*[1, 1]);

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

% clean up edges with morphological filters 
disk = strel('disk', morph_radius);
mask = imclose(imopen(mask, disk), disk);

% (optional) plot to facilitate parameter selection
if show
    figure()
    colormap(gray);
    subplot(3,2,1); imagesc(hue); title('hue'); set(gca,'XTick', [], 'YTick',[])
    subplot(3,2,2); imagesc(hue_mask); title('hue mask'); set(gca,'XTick', [], 'YTick',[])
    subplot(3,2,3); imagesc(value); title('value'); set(gca,'XTick', [], 'YTick',[])
    subplot(3,2,4); imagesc(value_mask); title('value mask'); set(gca,'XTick', [], 'YTick',[])
    subplot(3,2,5); imagesc(entropy); title('entropy'); set(gca,'XTick', [], 'YTick',[])
    subplot(3,2,6); imagesc(entropy_mask); title('entropy mask'); set(gca,'XTick', [], 'YTick',[])
    
    figure()    
    colormap(gray);
    subplot(2,1,1); imagesc(rgb2gray(rgb)); title('original'); set(gca,'XTick', [], 'YTick',[])
    subplot(2,1,2); imagesc(mask); title('mask'); set(gca,'XTick', [], 'YTick',[])    
end
