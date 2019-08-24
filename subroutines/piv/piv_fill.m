function [xx_fill, yy_fill, img_fill, mask_fill] = piv_fill(...
    xx, yy, img, mask, pad_r, pad_c)
% function [xx_fill, yy_fill, img_fill, mask_fill] = piv_fill(...
%   xx, yy, img, mask, pad_r, pad_c)
% 
% Fill non-sand regions of the input image by mirroring boundary pixels
%
% Arguments:
%   xx, yy: Vectors, coordinate axes for image columns, rows
%   img: 2D matrix, rectified/cropped grayscale image
%   mask: 2D matrix, logical mask of sand pixels in img
%   pad_r, pad_c: Scalar integers, number of pixels of padding to add in
%       row and column directions
% 
% Outputs:
%   xx_fill, yy_fill: Vectors, coordinate axes for filled image columns, rows
%   img_fill: 2D matrix, filled and padded grayscale image
%   mask_fill: 2D matrix, padded image mask matching (2D) size of img_fill
% %

% note: x-direction is padded with zeros and left masked. There is no
%   sane way to pad regions where sand grains enter or leave the frame.

% sanity checks
narginchk(6, 6);
validateattributes(xx, {'numeric'}, {'vector'});
validateattributes(yy, {'numeric'}, {'vector'});
validateattributes(img, {'double', 'single'}, {'2d', 'size', [numel(yy), numel(xx)]});
validateattributes(mask, {'logical'}, {'2d', 'size', [numel(yy), numel(xx)]});
validateattributes(pad_r, {'numeric'}, {'integer', 'nonnegative'});
validateattributes(pad_c, {'numeric'}, {'integer', 'nonnegative'});

fprintf('%s: fill and pad image and its mask by mirroring sand pixels\n', ...
    mfilename); 

% create padded arrays
img_fill = padarray(img, [pad_r, pad_c], NaN, 'both');  % FIXME: why tf am I nan padding here?
mask_fill = padarray(mask, [pad_r, pad_c], 0, 'both');  % FIXME: will want to return the mask afterwards, == has_sand?

% extend coordinate vectors
dy = yy(2) - yy(1);
yy_fill = [yy(1) - dy*(pad_r:-1:1), yy, yy(end) + dy*(1:pad_r)];

dx = xx(2) - xx(1);
xx_fill = [xx(1) - dx*(pad_c:-1:1), xx, xx(end) + dx*(1:pad_c)];
    
% find row index of the top and bottom boundaries
top_row = nan(size(xx_fill));
bot_row = nan(size(xx_fill));

for jj = 1:length(xx_fill)
    ii_bot = find(mask_fill(:, jj), 1, 'first');
    if ~isempty(ii_bot)
        bot_row(jj) = ii_bot;
    end
    ii_top = find(mask_fill(:, jj), 1, 'last');
    if ~isempty(ii_top)
        top_row(jj) = ii_top;
    end
end

% FIXME: what if, wait for it, this smooth line just was the mask all
%   along, I don't see a reason I could not store a pair of lines with 
%   non-integer mask edge in them during prep.

% smooth the upper boundary line
% note: lower boundary is *not* smooth due to the presence of the 
%   metal support in the image, leave it alone
smooth_num_pts = 10;
smooth_top_row = smooth(top_row, smooth_num_pts/length(top_row), 'lowess')';
smooth_top_row(isnan(top_row)) = nan;  % do not extrapolate! these regions have no pixels to mirror
top_row = smooth_top_row;  % rename and replace

% build index array by reflecting until all indices are within sand
% note: resulting coordinates are not integers, due to smoothing of the
%   upper boundary line
[nr, nc] = size(mask_fill);
bot_rows = repmat(bot_row, nr, 1);
top_rows = repmat(top_row, nr, 1);
[cols, rows] = meshgrid(1:nc, 1:nr);

has_sand = ... % only try to mirror where there is something to mirror
    ~isnan(top_rows) & ~isnan(bot_rows) & (top_rows - bot_rows > 3);
above_top = (rows > top_rows) & has_sand;

while any(above_top(:))

    % reflect at top boundary, may create inidices below bottom so update
    rows(above_top) = 2*top_rows(above_top) - rows(above_top);  % same as t-(r-t)
    below_bot = (rows < bot_rows) & has_sand;
    
    % reflect at bottom boundary, may create indices above top so update
    rows(below_bot) = 2*bot_rows(below_bot) - rows(below_bot);  % same as b+(b-r)
    above_top = (rows > top_rows) & has_sand;
    
end

% apply reflection by interpolating
% black out regions with no sand at all
% note: pixels within the mask are left unchanged, as desired
[jj, ii] = meshgrid(1:size(mask_fill, 2), 1:size(mask_fill, 1)); 
img_fill(:, :) = interp2(jj, ii, img_fill(:, :), cols, rows, 'bilinear');  % FIXME: higher order interpolants introduce NaNs

% black out regions with no sand at all
img_fill(~has_sand) = 0;

% sanity checks
assert(length(yy_fill) == size(img_fill, 1), 'Padded dimensions do not match padded image');
assert(length(xx_fill) == size(img_fill, 2), 'Padded dimensions do not match padded image');
