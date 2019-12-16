function [img_fill, mask_fill] = prep_fill_repeat(img, mask, skin_min, skin_max)
% function [img_fill, mask_fill] = prep_fill_repeat(img, mask, skin_min, skin_max)
%
% Fill non-sand regions of the input image by repeating a thin layer at sand boundary
%
% Arguments:
%   img: 2D matrix, rectified/cropped, and padded grayscale image
%   mask: 2D matrix, logical mask of sand pixels in img
%   skin_min: Optional integer, TODO
%   skin_sax: Optional integer, TODO
%
% Outputs:
%   img_fill: 2D matrix, filled grayscale image
%   mask_fill: 2D matrix, image mask matching size of img_fill
% %

% TODO: what if I ran this on the RGB panes, and equalized afterwards on a full image?

% TODO: consider what to do about left and right side padding

% TODO: consider what to do about regions of the image with no sand at all during buildup phase,
%   would be much better if they were filled so that I could drop (slow) masked cross correlation.

% TODO: consider special treatment for the center spacer, not clear what sort of padding will work
%   best in that region

narginchk(2, 4);
if nargin < 3; skin_min = 10; end
if nargin < 4; skin_max = 40; end

validateattributes(img, {'double', 'single'}, {'2d'});
validateattributes(mask, {'logical'}, {'2d', 'size', size(img)});
validateattributes(skin_min, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(skin_max, {'numeric'}, {'scalar', 'integer', 'positive', '>=', skin_min});


[nr, nc] = size(img);

% create repeatably random offset vector 

% TODO: replace with random permutations of possible skin widths, which should better assure we do
%   not repeat and get degenerate PIV results

% set seed, note that this modifies the global stream, which should be fine unless we require
% non-repeatable randomness somewhere downstream
rng(982198);

offsets  = nan(nr, 1);

start_idx = 1;
while start_idx < nr
    skin_depth = randi([skin_min, skin_max], 1);
    this = (skin_depth - 1):-1:1;
    end_idx = min(start_idx + length(this) - 1, length(offsets));
    offsets(start_idx:end_idx) = this(1:(end_idx - start_idx + 1));
    start_idx = end_idx + 1;
end
 
% find row index of the top and bottom boundaries
% TODO: rename to top_row_idx, etc, for clarity
top_row = nan(1, nc);
bot_row = nan(1, nc);

for jj = 1:nc
    ii_bot = find(mask(:, jj), 1, 'first');
    if ~isempty(ii_bot)
        bot_row(jj) = ii_bot;
    end
    ii_top = find(mask(:, jj), 1, 'last');
    if ~isempty(ii_top)
        top_row(jj) = ii_top;
    end
end

% TODO: does using a smoothed boundary really do anything? I am now skeptical, but it might be
%   useful just to force a (smooth) resampling in the filled region.

% smooth the upper boundary line
% note: lower boundary is *not* smooth due to the presence of the 
%   metal support in the image, leave it alone
% EXPERIMENT!
smooth_num_pts = 200;
smooth_top_row = smooth(top_row, smooth_num_pts/length(top_row), 'lowess')';
smooth_top_row(isnan(top_row)) = nan;  % do not extrapolate! these regions have no pixels to mirror
top_row = smooth_top_row;  % rename and replace

% EXPERIMENT: trim off the top layer
trim_pixels = 40;
top_row = max(1, top_row - trim_pixels);
 
% build index array
% note: resulting coordinates are not integers, due to smoothing of the upper boundary line
[cols, rows] = meshgrid(1:nc, 1:nr);

for jj = 1:nc
    
    if isnan(top_row(jj))
        % not enough to pad with, skip
        % TODO: handle case where there is something, but less than the skin depth
        continue
    end
    
    start_fill_idx = ceil(top_row(jj));
    num_fill = nr - start_fill_idx + 1;
    this_offsets = offsets + (start_fill_idx - top_row(jj));  % adjust so first offset mirrors boundary
    rows(start_fill_idx:end, jj) = top_row(jj) - this_offsets(1:num_fill);
    
end

for jj = 1:nc
    
    if isnan(bot_row(jj))
        % not enough to pad with, skip
        % TODO: handle case where there is something, but less than the skin depth
        continue
    end
    
    start_fill_idx = floor(bot_row(jj)); % TODO: this fails! we get a repeat value, try adding one then taking floor
    num_fill = start_fill_idx + 1;
    this_offsets = offsets + (bot_row(jj) - start_fill_idx);  % adjust so first offset mirrors boundary
    rows(1:num_fill, jj) = bot_row(jj) + this_offsets(num_fill:-1:1);
    
end

% TODO: assert something about spacing, always approx 1? always < 1? not sure what is right here.
 
% apply reflection by interpolating
% black out regions with no sand at all
% note: pixels within the mask are left unchanged, as desired
img_fill = img;
[jj, ii] = meshgrid(1:nc, 1:nr); 
% EXPERIMENT!
% img_fill(:, :) = interp2(jj, ii, img_fill(:, :), cols, rows, 'bilinear');  % FIXME: higher order interpolants introduce NaNs
img_fill(:, :) = interp2(jj, ii, img_fill(:, :), cols, rows, 'cubic');  % FIXME: higher order interpolants introduce NaNs
mask_fill = ~isnan(img_fill);