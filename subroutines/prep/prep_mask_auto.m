function clean_mask = prep_mask_auto(...
    rgb, hue_lim, value_lim, entropy_lim, entropy_len, view, show)
% function clean_mask = prep_mask_auto(hsv, hue_lim, value_lim, entropy_lim, entropy_len, ...
%                     morph_open_rad, morph_erode_rad, view, show)
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
%   view =  string, specify 'side' or 'top' image view
%
%   show = Scalar, logical, set to 1 (true) to plot the mask bands, used to
%       facilitate the parameter selection process, default = false.
% 
%   mask = 2D matrix, logical, true where there is sand and false
%       elsewhere.
% % 

% set default values
if nargin < 7; show = false; end

% check for sane arguments, set default values
narginchk(6, 7);
validateattributes(rgb, {'numeric'}, {'3d'});
validateattributes(hue_lim, {'double'}, {'vector', 'numel', 2, '>=', 0, '<=', 1});
validateattributes(value_lim, {'double'}, {'vector', 'numel', 2, '>=', 0, '<=', 1});
validateattributes(entropy_lim, {'double'}, {'vector', 'numel', 2, '>=', 0, '<=', 1});
validateattributes(entropy_len, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(view, {'char'}, {});
validateattributes(show, {'numeric', 'logical'}, {'scalar'});

% get size constants
nr = size(rgb, 1);
nc = size(rgb, 2);

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
ridx = 1:nr;
cidx = 1:nc;
dr = [1, 1, 0, 0];
dc = [1, 0, 0, 1];
for ii = 1:4
    
    wall = true([nr, nc]+1);
    wall(ridx+dr(ii), cidx+dc(ii)) = raw_mask;
    wall = imfill(wall, 'holes');
    raw_mask = wall(ridx+dr(ii), cidx+dc(ii));
end
 
% extract largest connected object
object_label = bwlabel(raw_mask);
largest_object = mode(object_label(object_label>0));
raw_mask = object_label == largest_object;

% create a mask without holes from the first and last pixels per column
rough_mask = false([nr, nc]);
for jj = 1:nc
    column = raw_mask(:, jj);
    ii_min = find(column, 1, 'first');
    ii_max = find(column, 1, 'last');
    rough_mask(ii_min:ii_max, jj) = true;
end

% clean up mask edges
threshold_percentile =  25; % note: this threshold percentile works so far

if strcmp(view, 'side')
    % drop pixels from upper sand boundary based on a threshold brightness
    %   (value) computed from the rough mask, this effectively removes the
    %   "halo" of background pixels in the rough mask
    
    % check assumption that wedge top is at high indices and bottom at low
    assert(sum(rough_mask(1:10, :), 'all') > sum(rough_mask(end-10:end, :), 'all'), ...
        'found more sand at high indices than at low, but expected wedge top at high indices');
        
    % get raw threshold value in each column
    % note: smooth threshold values to get a stable estimate
    data_idx = find(max(rough_mask));  % find columns that have some sand
    threshold_value = zeros(1, nc);
    min_idx = ones(1, nc);
    max_idx = ones(1, nc);
    for jj = data_idx  % leave min/max idx as 1 for columns with no sand 
        min_idx(jj) = find(rough_mask(:, jj), 1, 'first');
        max_idx(jj) = find(rough_mask(:, jj), 1, 'last');
        threshold_value(jj) = prctile(value(min_idx(jj):max_idx(jj), jj), threshold_percentile);
    end
    threshold_value = smooth(threshold_value, 0.3, 'lowess');
    
    % find the new top boundary 
    masked_value = value;
    masked_value(~rough_mask) = NaN;
    data_idx =  find(max(masked_value > repmat(threshold_value', size(masked_value, 1), 1)));  % find columns that have pixels above the threshold
    threshold_max_idx = ones(1, nc);
    for jj = data_idx
        threshold_max_idx(jj) = find(masked_value(:, jj) >= threshold_value(jj), 1, 'last');
    end

    % replace outliers with interpolated estimate
    outliers = false(size(threshold_max_idx));
    outliers(data_idx) = isoutlier(...
        max_idx(data_idx) - threshold_max_idx(data_idx), 'gesd')';  %  generalized extreme Studentized deviate test 
    cols = 1:nc;
    threshold_max_idx(outliers) = interp1(...
        cols(~outliers)', threshold_max_idx(~outliers)', ...
        cols(outliers), 'pchip', 0);  % extrapolated values set to 0 (below bottom)
    
    % apply a touch of smoothing to the new boundary
    threshold_max_idx = round(smooth(threshold_max_idx, 5/nc, 'lowess'));
    
    % drop edge pixels below threshold values
    clean_mask = false([nr, nc]);
    for jj = 1:nc
        clean_mask(min_idx(jj):threshold_max_idx(jj), jj) = true;
    end
    
elseif strcmp(view, 'top')
    % not clear how top view should be handled, warn and skip for now
    warning('edge cleanup for top view is not implemented');
    clean_mask = rough_mask;
end

% report percentage masked
pct_sand = 100*sum(clean_mask(:))/numel(clean_mask);
fprintf('%s: %.0f%% sand, %.0f%% background\n', mfilename, pct_sand, 100-pct_sand);

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
    subplot(2,1,2); imagesc(clean_mask); title('mask'); set(gca,'XTick', [], 'YTick',[])
    
    figure()
    rough_mask_bnd = bwboundaries(rough_mask);
    rough_mask_x = rough_mask_bnd{1}(:,2);
    rough_mask_y = rough_mask_bnd{1}(:,1);
    clean_mask_bnd = bwboundaries(clean_mask);
    clean_mask_x = clean_mask_bnd{1}(:,2);
    clean_mask_y = clean_mask_bnd{1}(:,1);
    imagesc(value); colormap('gray'); hold on;
    plot(rough_mask_x, rough_mask_y, '-r');
    plot(clean_mask_x, clean_mask_y, '-b');
    axis equal tight
    set(gca, 'YDir', 'normal', 'XTick', [], 'YTick',[]);
    title('Clean (blue) and rough (red) boundary lines');
end
