function mask = yalebox_prep_mask_auto(hsv, hue_lim, val_lim, entr_lim, med_win, entr_win, morph_rad, show)
% function mask = yalebox_prep_mask_auto(hsv, hue_lim, val_lim, entr_lim, med_win, entr_win, morph_rad, show)
%
% Create a logical mask for a color image that is TRUE where there is sand and
% FALSE elsewhere. This can be used to remove (set to 0) the background in a
% image prior to PIV analysis or other applications. Sand is identified by
% remapping to HSV colorspace, filtering and thresholding the "hue" and "value"
% bands.
%
% Arguments:
%
%   hsv = 3D matrix, double, Color image in HSV colorspace.
%
%   hue_lim = 2-element vector, double, range [0, 1]. [minimum, maximum] HSV
%     "hue" included as sand in the mask.
%
%   val_lim = 2-element vector, double, range [0,1]. [minimum, maximum] HSV
%     "value" included as sand in the mask.
%
%   entr_lim = 2-element vector, double, range [0, 1]. [minimum, maximum]
%     entropy included as sand in the mask. 
% 
%   med_win = Scalar, integer, window size in pixels for median filter.
%
%   entr_win = scalar, integer, window size in pixels for entropy filter.
%
%   morph_rad = scalar, double, radius of disk structuring element used in
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
assert(isa(hsv, 'double') && size(hsv,3) == 3, ...
    'hsv is not a HSV image');
assert(numel(hue_lim)==2 && min(size(hue_lim)) == 1, ...
    'hue_lim is not a 2-element vector');
assert(min(hue_lim) >= 0 && max(hue_lim) <= 1, ...
    'hue_lim is not in the range [0,1]');
assert(numel(val_lim)==2 && min(size(val_lim)) == 1, ...
    'val_lim is not a 2-element vector');
assert(min(val_lim) >= 0 && max(val_lim) <= 1, ...
    'val_lim is not in the range [0,1]');
assert(numel(entr_lim)==2 && min(size(entr_lim)) == 1, ...
    'entr_lim is not a 2-element vector');
assert(min(entr_lim) >= 0 && max(entr_lim) <= 1, ...
    'entr_lim is not in the range [0,1]');
assert(numel(med_win) == 1 && med_win > 0 && ...
    round(med_win) == med_win, ...
    'med_win is not a positive integer');
assert(numel(entr_win) == 1 && entr_win > 0 && ...
    round(entr_win) == entr_win, ...
    'entr_win is not a positive integer');
assert(numel(morph_rad) == 1 && morph_rad > 0, ...
    'morph_rad is not a positive scalar');
if nargin == 7; show = false; end
assert(numel(show) == 1 && ismember(show, [0, 1]), ...
    'show is not 1 or 0');

% get hue, value, and entropy, normalized to the range [0, 1]
hue = hsv(:,:,1);
value = hsv(:,:,3);
entropy = entropyfilt(value, true(entr_win));

hue = hue-min(hue(:)); hue = hue./max(hue(:));
value = value-min(value(:)); value = value./max(value(:));
entropy = entropy-min(entropy(:)); entropy = entropy./max(entropy(:));

% threshold bands
hue_mask = hue >= hue_lim(1) & hue <= hue_lim(2);
value_mask = value >= val_lim(1) & value <= val_lim(2);
entropy_mask = entropy >= entr_lim(1) & entropy <= entr_lim(2);

% create mask
mask = hue_mask & value_mask & entropy_mask;

% smooth, preseving edges with median filter
mask = medfilt2(mask, med_win*[1, 1]);

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
disk = strel('disk', morph_rad);
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
    subplot(2,1,1); imagesc(hsv(:,:,3)); title('original'); set(gca,'XTick', [], 'YTick',[])
    subplot(2,1,2); imagesc(mask); title('mask'); set(gca,'XTick', [], 'YTick',[])    
end
