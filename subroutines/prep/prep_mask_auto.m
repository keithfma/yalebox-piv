function mask = prep_mask_auto(rgb, hue_lim, value_lim, entropy_lim, entropy_len, ...
                    morph_open_rad, morph_erode_rad, show)
% function mask = prep_mask_auto(hsv, hue_lim, value_lim, entropy_lim, entropy_len, ...
%                     morph_open_rad, morph_erode_rad, show)
%
% Create a logical mask for a color image that is TRUE where there is sand and
% FALSE elsewhere. This can be used to remove (set to 0) the background in a
% image prior to PIV analysis or other applications. Sand is identified by
% thresholding "hue", "value" and "entropy" bands.
%
% Arguments:
%
%   rgb = 3D matrix, color image in RGB colorspace.
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
%   entropy_len = scalar, integer, window size in pixels for entropy filter.
%
%   morph_open_rad = scalar, double, radius of disk structuring element used in
%     mophological opening filter. 
%
%   morph_erode_rad = scalar, double, radius of disk structuring element used in
%     mophological erosion filter.
%
%   show = Scalar, logical, set to 1 (true) to plot the mask bands, used to
%       facilitate the parameter selection process, default = false.
% 
%   mask = 2D matrix, logical, true where there is sand and false
%       elsewhere.
% % 

% set default values
if nargin < 8; show = false; end

% check for sane arguments, set default values
narginchk(7, 8);
validateattributes(rgb, {'numeric'}, {'3d'});
validateattributes(hue_lim, {'double'}, {'vector', 'numel', 2, '>=', 0, '<=', 1});
validateattributes(value_lim, {'double'}, {'vector', 'numel', 2, '>=', 0, '<=', 1});
validateattributes(entropy_lim, {'double'}, {'vector', 'numel', 2, '>=', 0, '<=', 1});
validateattributes(entropy_len, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(morph_open_rad, {'numeric'}, {'scalar', '>=', 0});
validateattributes(morph_erode_rad, {'numeric'}, {'scalar', '>=', 0});
validateattributes(show, {'numeric', 'logical'}, {'scalar'});

% get hue, value, and entropy, normalized to the range [0, 1]
hsv = rgb2hsv(rgb);
hue = hsv(:,:,1);
value = hsv(:,:,3);
entropy = entropyfilt(value, true(entropy_len));

hue = hue-min(hue(:)); hue = hue./max(hue(:));
value = value-min(value(:)); value = value./max(value(:));
entropy = entropy-min(entropy(:)); entropy = entropy./max(entropy(:));

% threshold bands
hue_mask = hue >= hue_lim(1) & hue <= hue_lim(2);
value_mask = value >= value_lim(1) & value <= value_lim(2);
entropy_mask = entropy >= entropy_lim(1) & entropy <= entropy_lim(2);

% create mask
mask = hue_mask & value_mask & entropy_mask;

% fill holes along edges (wall off one corner, fill, repeat)
ridx = 1:size(mask, 1);
cidx = 1:size(mask, 2);
dr = [1, 1, 0, 0];
dc = [1, 0, 0, 1];
for ii = 1:4
    wall = true(size(mask)+1);
    wall(ridx+dr(ii), cidx+dc(ii)) = mask;
    wall = imfill(wall, 'holes');
    mask = wall(ridx+dr(ii), cidx+dc(ii));
end

% extract largest connected object
object_label = bwlabel(mask);
largest_object = mode(object_label(object_label>0));
mask = object_label == largest_object;

% DEBUG: try "singed" mask to improve reliability

% clean up edges with morphological filters
% ...remove noise (open), then trim boundary (erode)
if morph_open_rad > 0
    mask = imopen(mask, strel('disk', morph_open_rad));
end

threshold = 0.5;

hsv = rgb2hsv(rgb);
val = hsv(:, :, 3);

for cc = 1:size(mask, 2)
    not_nan = ~isnan(mask(:, cc));
    
    
    median_val = median(val(not_nan, cc));
    above_threshold = val(:, cc) > median_val*threshold;    
    
    keep = not_nan & above_threshold;
    
    min_idx = find(keep, 1, 'first');
    max_idx = find(keep, 1, 'last');
    mask(1:(min_idx-1), cc) = false;
    mask((max_idx+1):end, cc) = false;
end

% if morph_erode_rad > 0
%     mask = imerode(mask, strel('disk', morph_erode_rad));
% end

% END DEBUG

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

% (optional) report percentage masked
pct_sand = 100*sum(mask(:))/numel(mask);
fprintf('%s: %.0f%% sand, %.0f%% background\n', mfilename, pct_sand, 100-pct_sand);
