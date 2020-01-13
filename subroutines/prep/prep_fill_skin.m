function [img_fill, mask_fill] = prep_fill_skin(img, mask, skin_min, skin_max)
% function img_fill = prep_fill(img, mask)
% 
% Fill non-sand regions of the input image by mirroring across sand boundary
%
% Arguments:
%   img: 2D matrix, rectified/cropped, and padded grayscale image
%   mask: 2D matrix, logical mask of sand pixels in img
%   skin_min: Optional integer, minimum thickness of skin layer to mirror when filling, 
%       will choose a random series of thicknesses between 'skin_min' and 'skin_max'
%   skin_sax: Optional integer, maximum thickness of skin layer to mirror when filling,
%       will choose a random series of thicknesses between 'skin_min' and 'skin_max'
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

narginchk(2, 4);
if nargin < 3; skin_min = 20; end
if nargin < 4; skin_max = 30; end

validateattributes(img, {'double', 'single'}, {'2d'});
validateattributes(mask, {'logical'}, {'2d', 'size', size(img)});
validateattributes(skin_min, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(skin_max, {'numeric'}, {'scalar', 'integer', 'positive', '>=', skin_min});


[nr, nc] = size(img);

% create repeatably random offset vector 

% set seed, note that this modifies the global stream, which should be fine unless we require
% non-repeatable randomness somewhere downstream
rng(982198);

skin_depth_idx = 1;
skin_depths = []; 

offsets  = nan(nr, 1);
offsets_dist = 0:(nr-1);  % first element of the offset vector corresponds to the boundary

start_idx = 1;
while start_idx < nr
    
    if skin_depth_idx > length(skin_depths)
        skin_depths =  randperm(skin_max - skin_min + 1) + skin_min -1;  % random permutation of ints in range [skin_min, skin_max], inclusive
        skin_depth_idx = 1;
    end
    skin_depth = skin_depths(skin_depth_idx);
    skin_depth_idx = skin_depth_idx + 1;
        
    this = [0:skin_depth, (skin_depth - 1):-1:1];
    end_idx = min(start_idx + length(this) - 1, length(offsets));
    offsets(start_idx:end_idx) = this(1:(end_idx - start_idx + 1));
    start_idx = end_idx + 1;
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

% smooth the upper boundary line
% note: lower boundary is *not* smooth due to the presence of the 
%   metal support in the image, leave it alone
smooth_num_pts = 10;
smooth_top_row_idx = smooth(top_row_idx, smooth_num_pts/length(top_row_idx), 'lowess')';
smooth_top_row_idx(isnan(top_row_idx)) = nan;  % do not extrapolate! these regions have no pixels to mirror
top_row_idx = smooth_top_row_idx;  % rename and replace
 
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

% % DEBUG 1111: review the index array to search for hard boundaries
% figure
% imagesc(rows);
% caxis([265, 285]);
% colorbar
% hold on;
% plot(top_row_idx, 'Marker', '.');
% axis equal
% set(gca, 'XLim', [3688, 3816], 'YLim', [266, 332]);
% keyboard
% % /DEBUG 1111

% apply reflection by interpolating
% black out regions with no sand at all
% note: pixels within the mask are left unchanged, as desired
img_fill = img;
[jj, ii] = meshgrid(1:nc, 1:nr); 
img_fill(:, :) = interp2(jj, ii, img_fill(:, :), cols, rows, 'bilinear');  % FIXME: higher order interpolants introduce NaNs
mask_fill = ~isnan(img_fill);
