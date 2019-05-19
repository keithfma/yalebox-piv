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
raw_mask = hue_mask & value_mask & entropy_mask;

% fill holes along edges (wall off one corner, fill, repeat)
ridx = 1:size(raw_mask, 1);
cidx = 1:size(raw_mask, 2);
dr = [1, 1, 0, 0];
dc = [1, 0, 0, 1];
for ii = 1:4
    wall = true(size(raw_mask)+1);
    wall(ridx+dr(ii), cidx+dc(ii)) = raw_mask;
    wall = imfill(wall, 'holes');
    raw_mask = wall(ridx+dr(ii), cidx+dc(ii));
end
 
% extract largest connected object
object_label = bwlabel(raw_mask);
largest_object = mode(object_label(object_label>0));
raw_mask = object_label == largest_object;

% create a mask without holes from the first and last pixels per column
rough_mask = false(size(raw_mask));
for jj = 1:size(raw_mask, 2)
    column = raw_mask(:, jj);
    ii_min = find(column, 1, 'first');
    ii_max = find(column, 1, 'last');
    rough_mask(ii_min:ii_max, jj) = true;
end

% clean up mask edges
% note: low threshold percentile works, avoid parameterizing if possible
mask = singe(value, rough_mask, 25);

% TODO; singing working well at upper bnd but not lower, suspect a
%   different threshold is needed. One approach would be to do stats on 
%   a layer of fixed width from the boundary only, and run the top and
%   bottom separately.

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
    
    figure()
    rough_mask_bnd = bwboundaries(rough_mask);
    rough_mask_x = rough_mask_bnd{1}(:,2);
    rough_mask_y = rough_mask_bnd{1}(:,1);
    mask_bnd = bwboundaries(mask);
    mask_x = mask_bnd{1}(:,2);
    mask_y = mask_bnd{1}(:,1);
    imagesc(value); colormap('gray'); hold on;
    plot(rough_mask_x, rough_mask_y, '-r');
    plot(mask_x, mask_y, '-b');
    axis equal tight
    set(gca, 'YDir', 'normal', 'XTick', [], 'YTick',[]);
    title('Final (blue) and rough (red) boundary lines');
end

% (optional) report percentage masked
pct_sand = 100*sum(mask(:))/numel(mask);
fprintf('%s: %.0f%% sand, %.0f%% background\n', mfilename, pct_sand, 100-pct_sand);


function new_mask = singe(value, mask, threshold_percentile)
% function mask = singe(rgb, mask, threshold_percentile)
%
% Improve mask edge quality by burning off ("singing") low-brightness pixels
%
% Arguments:
%   value: 
%   mask: boolean matrix, true for sand pixels, false elsewhere
%   threshold percentile: float, range 0-100, threshold below which pixels
%       are "singed" out of the input mask
%
% Returns:
%   mask: singed version of the input mask
% %

% TODO: even with stable threshold estimate, there are some columns that 
%   are deeply singed next to other that are not. This yields an irregular
%   comb shape on mask edges, especially the bottom. Try smoothing the
%   singed boundary, or something else if that fails

% get raw threshold value for each column
threshold = zeros(1, size(value, 2));
for jj = 1:size(mask, 2)
    % compute first and last index of mask in each column
    min_idx = find(mask(:, jj), 1, 'first');
    max_idx = find(mask(:, jj), 1, 'last');
    % get raw threshold value
    % note: check assumption that mask is all true between min and max
    assert(all(mask(min_idx:max_idx, jj)), 'mask has holes in this column, bad');
    threshold(jj) = prctile(value(min_idx:max_idx, jj), threshold_percentile);
end

% smooth threshold value to get a stable estimate
threshold = smooth(threshold, 0.3, 'loess');

% singe edge pixels below threshold
masked_val = value;
masked_val(~mask) = NaN;
new_mask = false(size(mask));

for jj = 1:size(mask, 2)
    min_idx = find(masked_val(:, jj) >= threshold(jj), 1, 'first');
    max_idx = find(masked_val(:, jj) >= threshold(jj), 1, 'last');
    new_mask(min_idx:max_idx, jj) = true;
end


