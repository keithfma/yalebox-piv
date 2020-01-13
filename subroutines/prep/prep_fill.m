function [img_fill, mask_fill] = prep_fill_skin(img, mask, skin_min, skin_max, bnd_smooth_window, mirror)
% function [img_fill, mask_fill] = prep_fill_skin(img, mask, skin_min, skin_max, bnd_smooth_window, mirror)
% 
% Fill non-sand regions of the input image by mirroring across sand boundary
%
% Arguments:
%   img: 2D matrix, rectified/cropped, and padded grayscale image
%   mask: 2D matrix, logical mask of sand pixels in img
%   skin_min: Optional integer, minimum thickness of skin layer to use when filling, 
%       will choose a random series of thicknesses between 'skin_min' and 'skin_max'
%   skin_max: Optional integer, maximum thickness of skin layer to use when filling,
%       will choose a random series of thicknesses between 'skin_min' and 'skin_max'
%   bnd_smooth_window: width of the (loess) smoothing window to apply to the upper 
%       boundary before computing the skin layer, set to 0 to disable smoothing
%   mirror: Optional logical, set True to mirror the skin layer when filling, or
%       false (default) to repeat it
%   show: Logical, set true to display results for visual debugging, default is false
%
% Outputs:
%   img_fill: 2D matrix, filled grayscale image
%   mask_fill: 2D matrix, image mask matching size of img_fill
% %

% TODO: consider what to do about left and right side padding

% TODO: consider what to do about regions of the image with no sand at all during buildup phase,
%   would be much better if they were filled so that I could drop (slow) masked cross correlation.

% TODO: consider special treatment for the center spacer, not clear what sort of padding will work
%   best in that region

% TODO 1111: make skin_min, skin_max, and bnd_smooth_window required arguments
narginchk(2, 7);
if nargin < 3; skin_min = 3; end
if nargin < 4; skin_max = 10; end
if nargin < 5; bnd_smooth_window = 10; end
if nargin < 6; mirror = false; end
if nargin < 7; show = false; end

validateattributes(img, {'double', 'single'}, {'2d'});
validateattributes(mask, {'logical'}, {'2d', 'size', size(img)});
validateattributes(skin_min, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(skin_max, {'numeric'}, {'scalar', 'integer', 'positive', '>=', skin_min});
validateattributes(bnd_smooth_window, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
validateattributes(mirror, {'logical'}, {'scalar'});

fprintf('%s: fill image: skin=[%d, %d], bnd_smooth_window=%d, mirror=%d\n', ...
    mfilename, skin_min, skin_max, bnd_smooth_window, mirror); 

[nr, nc] = size(img);

% set seed, note that this modifies the global stream, which should be fine unless we require
% non-repeatable randomness somewhere downstream
rng(982198);

% create repeatably random offset vector 
if mirror
    [offsets, offsets_dist] = mirror_offsets(nr, skin_min, skin_max);
else:
    [offsets, offsets_dist] = repeat_offsets(nr, skin_min, skin_max);
end

% find row index of the top and bottom boundaries
top_row_idx = nan(1, nc);
bot_row_idx = nan(1, nc);

for jj = 1:nc
    ii_bot = find(mask(:, jj), 1, 'first');
    if ~isempty(ii_bot)
        bot_row_idx(jj) = ii_bot;
    end
    ii_top = find(mask(:, jj), 1, 'last');
    if ~isempty(ii_top)
        top_row_idx(jj) = ii_top;
    end
end

if bnd_smooth_window > 0
    % smooth the upper boundary line
    % note: lower boundary is *not* smooth due to the presence of the 
    %   metal support in the image, leave it alone
    smooth_top_row_idx = smooth(top_row_idx, bnd_smooth_window/length(top_row_idx), 'lowess')';
    smooth_top_row_idx(isnan(top_row_idx)) = nan;  % do not extrapolate! these regions have no pixels to mirror
    top_row_idx = smooth_top_row_idx;  % rename and replace
end
 
% build index array
% note: resulting coordinates are not integers, due to smoothing of the upper boundary line
[cols, rows] = meshgrid(1:nc, 1:nr);

for jj = 1:nc
    
    if isnan(top_row_idx(jj))
        % not enough to pad with, skip
        % TODO: handle case where there is something, but less than the skin depth
        continue
    end
    
    to_fill = rows(:, jj) > top_row_idx(jj);
    num_fill = sum(to_fill);
    
    % get offsets by interpolating into the offsets vector based on distance to boundary
    dist = rows(to_fill, jj) - top_row_idx(jj);
    this_offsets = interp1(offsets_dist(1:num_fill), offsets(1:num_fill), dist);

    rows(to_fill, jj) = top_row_idx(jj) - this_offsets;

end

% TODO 1111: use same method as top fill here, abstract away to a helper

for jj = 1:nc
    
    if isnan(bot_row_idx(jj))
        % not enough to pad with, skip
        % TODO: handle case where there is something, but less than the skin depth
        continue
    end
    
    start_fill_idx = floor(bot_row_idx(jj)); % TODO: this fails! we get a repeat value, try adding one then taking floor
    num_fill = start_fill_idx + 1;
    this_offsets = offsets + (bot_row_idx(jj) - start_fill_idx);  % adjust so first offset mirrors boundary
    rows(1:num_fill, jj) = bot_row_idx(jj) + this_offsets(num_fill:-1:1);
    
end

% apply reflection by interpolating
% black out regions with no sand at all
% note: pixels within the mask are left unchanged, as desired
img_fill = img;
[jj, ii] = meshgrid(1:nc, 1:nr); 
img_fill(:, :) = interp2(jj, ii, img_fill(:, :), cols, rows, 'bilinear');  % FIXME: higher order interpolants introduce NaNs
mask_fill = ~isnan(img_fill);

if show
    hf = figure;
    hf.Units = 'normalized';
    hf.Position = [0, hf.Position(2), 1, hf.Position(4)];
    
    imagesc(img_fill);
    hold on
    bnds = bwboundaries(mask);  % not fill, original mask
    for ii = 1:length(bnds)
        bnd = bnds{ii};
        plot(bnd(:, 2), bnd(:, 1), 'r', 'LineWidth', 2);
    end
            
    colormap('gray');
    hax = gca;
    hax.YDir = 'normal';
    hax.YTick = [];
    hax.XTick = [];
    axis equal tight
    title(sprintf('Image filled with method "%s"', method'));
end

end


function [offs, offsets_dist] = mirror_offsets(num_rows, min_depth, max_depth)
% function [offsets, offsets_dist] = mirror_offsets(num_rows, min_depth, max_depth)
% 
% Create vector of offsets (from boundary inwards) for filled region by mirroring the skin layer
%
% Arguments:
%   num_rows: number of rows in the image to be padded
%   min_depth: see prep_fill_skin argument skin_min
%   max_depth: see prep_fill_skin argument skin max
% 
% Outputs:
%   offs = vector of offsets from the boundary for the padded region (integer)
%   offs_dist = vector of distances from the boundary corresponding to 'offs' vector 
% % 

depth_idx = 1;
depths = []; 

offs  = nan(num_rows, 1);
offs_dist = 0:(num_rows-1);  % first element of the offset vector corresponds to the boundary

start_idx = 1;
while start_idx < num_rows
    
    if depth_idx > length(depths)
        depths =  randperm(max_depth - min_depth + 1) + min_depth -1;  % random permutation of ints in range [min_depth, max_depth], inclusive
        depth_idx = 1;
    end
    skin_depth = depths(depth_idx);
    depth_idx = depth_idx + 1;
        
    this = [0:skin_depth, (skin_depth - 1):-1:1];
    end_idx = min(start_idx + length(this) - 1, length(offs));
    offs(start_idx:end_idx) = this(1:(end_idx - start_idx + 1));
    start_idx = end_idx + 1;
end

end

 
function [offs, offsets_dist] = repeat_offsets(num_rows, min_depth, max_depth)
% function [offsets, offsets_dist] = repeat_offsets(num_rows, min_depth, max_depth)
% 
% Create vector of offsets (from boundary inwards) for filled region by repeating the skin layer
%
% Arguments:
%   num_rows: number of rows in the image to be padded
%   min_depth: see prep_fill_skin argument skin_min
%   max_depth: see prep_fill_skin argument skin max
% 
% Outputs:
%   offs = vector of offsets from the boundary for the padded region (integer)
%   offs_dist = vector of distances from the boundary corresponding to 'offs' vector 
% % 

offs_parts = cell(0);
offs_dist_parts = cell(0);

max_dist = 0;
while max_dist < nr
    depth = randi([min_depth, max_depth], 1);
    
    this_offs = (depth - 1):-1:0;
    
    this_offs_dist = max_dist + (0:length(this_offs) - 1);
    this_offs_dist(1) = this_offs_dist(1) + eps(this_offs_dist(1));
    
    offs_parts{end + 1} = this_offs;  %#ok!
    offs_dist_parts{end + 1} = this_offs_dist;  %#ok!
    
    max_dist = this_offs_dist(end);
end

offs = cell2mat(offs_parts);
offs_dist = cell2mat(offs_dist_parts);

end

