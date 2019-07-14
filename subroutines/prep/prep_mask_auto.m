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

% TODO: could I write a tool that lets you click in one boundary, then
%   selects the mask parameters to best match it? Bet I could...in fact, 
%   this sounds a lot like the label/decision tree model I build a while
%   back

% TODO: consider order of operations to make sure mirror pad still works a
%   trick on this - singing may no longer be needed

% TODO: consider splitting the prep step into rectify, mask, and perhaps
%   mirror steps, doing so might allow me to get rid of the script. A
%   function to update a parameter value would help here too and make the
%   workflow more functional.

% debug: turn in-plane masking on or off for comparison
% do_in_plane_mask = true;
do_in_plane_mask = false;

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

pct_sand = 100*sum(rough_mask(:))/numel(rough_mask);
fprintf('%s: rough mask: %.1f%% sand\n', mfilename, pct_sand);

% remove out-of-plane pixels from top of wedge in side view
% note: these pixels manifest as a bright upper layer where the lighting
%   shines on the top of the sand
if do_in_plane_mask
    if strcmp(view, 'side')

        % constants
        layer_thickness = 70; % FIXME: this will fail when the top boundary is too low...
        smooth_width = 300; % FIXME: should probably be an input argument
        smooth_height = 1;
        threshold_value = 0.05;  % FIXME: select this automatically

        % check assumption that wedge top is at high indices and bottom at low
        assert(sum(rough_mask(1:10, :), 'all') > sum(rough_mask(end-10:end, :), 'all'), ...
            'found more sand at high indices than at low, but expected wedge top at high indices');

        % copy out the upper sand layer to a square array
        value_idx = reshape(1:numel(value), size(value));
        layer_idx = nan(layer_thickness, size(value, 2));
        for jj = 1:size(value, 2)
            top = find(rough_mask(:, jj), 1, 'last');
            layer_idx(:, jj) = value_idx((top - layer_thickness+1):top, jj);
        end
        layer = reshape(value(layer_idx), size(layer_idx)); 

        % detrend x-direction lighting gradients in the layer array
        % FIXME: may need to rescale to some known range, unless automatic
        %   threshold selection makes this irrelevant
        trend = smooth(median(layer, 1), 0.2, 'lowess');
        layer = layer - repmat(trend', [size(layer, 1), 1]);

        % pad, smooth, and unpad
        kernel = fspecial('average', [smooth_height, smooth_width]);
        padded_layer = padarray(layer, [smooth_height, smooth_width], 'symmetric', 'both');
        padded_layer = imfilter(padded_layer, kernel);
        layer = padded_layer((1+smooth_height):(end-smooth_height), (1+smooth_width):(end-smooth_width));

        % create mask from layer by filling upwards
        layer_mask = true(size(layer));
        for jj = 1:size(layer, 2)
            bottom = find(layer(:, jj) >= threshold_value, 1, 'first');
            layer_mask((bottom+1):end, jj) = false;
        end

        % expand mask from layer size to full image size
        layer_mask_full = true(size(value));
        layer_mask_full(layer_idx) = layer_mask;

        % combine with rough mask
        in_plane_mask = rough_mask & layer_mask_full; 

    elseif strcmp(view, 'top') || ~do_in_plane_mask
        warning('no out-of-plane masking for top view');
        in_plane_mask = rough_mask;
    end

else
    warning('out-of-plane masking disabled');
    in_plane_mask = rough_mask;
end    

pct_sand = 100*sum(in_plane_mask(:))/numel(in_plane_mask);
fprintf('%s: layer mask: %.1f%% sand\n', mfilename, pct_sand);

% clean up mask edges
% note: thresholding commonly leaves a "halo" of non-sand pixels within
%   the mask, which have been problematic for image padding routines

if strcmp(view, 'side')
    % drop pixels from upper sand boundary based on a threshold brightness
    %   (value) computed from the layer mask, this effectively removes the
    %   "halo" of background pixels in the layer mask

    threshold_percentile =  25; % note: this threshold percentile works so far
    
    % check assumption that wedge top is at high indices and bottom at low
    assert(sum(in_plane_mask(1:10, :), 'all') > sum(in_plane_mask(end-10:end, :), 'all'), ...
        'found more sand at high indices than at low, but expected wedge top at high indices');
        
    % get raw threshold value in each column
    % note: smooth threshold values to get a stable estimate
    threshold_value = zeros(1, nc);
    for jj = 1:nc
        min_idx = find(in_plane_mask(:, jj), 1, 'first');
        max_idx = find(in_plane_mask(:, jj), 1, 'last');
        threshold_value(jj) = prctile(value(min_idx:max_idx, jj), threshold_percentile);
    end
    threshold_value = smooth(threshold_value, 0.3, 'loess');
    
    % drop edge pixels below threshold values
    masked_value = value; masked_value(~in_plane_mask) = NaN;
    clean_mask = false([nr, nc]);
    for jj = 1:nc
        min_idx = find(in_plane_mask(:, jj), 1, 'first');  % note: same as above, lower edge already clean cleaned
        max_idx = find(masked_value(:, jj) >= threshold_value(jj), 1, 'last');
        clean_mask(min_idx:max_idx, jj) = true;
    end
    
elseif strcmp(view, 'top')
    % not clear how top view should be handled, warn and skip for now
    warning('edge cleanup for top view is not implemented');
    clean_mask = in_plane_mask;

end
pct_sand = 100*sum(clean_mask(:))/numel(clean_mask);
fprintf('%s: clean mask: %.1f%% sand\n', mfilename, pct_sand);

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

    in_plane_mask_bnd = bwboundaries(in_plane_mask);
    in_plane_mask_x = in_plane_mask_bnd{1}(:,2);
    in_plane_mask_y = in_plane_mask_bnd{1}(:,1);
    
    clean_mask_bnd = bwboundaries(clean_mask);
    clean_mask_x = clean_mask_bnd{1}(:,2);
    clean_mask_y = clean_mask_bnd{1}(:,1);
    
    imagesc(value); colormap('gray'); hold on;
    plot(rough_mask_x, rough_mask_y, '-r');
    plot(in_plane_mask_x, in_plane_mask_y, '-b');
    plot(clean_mask_x, clean_mask_y, '-g');
    axis equal tight
    set(gca, 'YDir', 'normal', 'XTick', [], 'YTick',[]);
    title('Clean (green), in-plane (blue), and rough (red) boundary lines');
end
