function [xx_fill, yy_fill, img_fill, mask_fill] = piv_fill(...
    xx, yy, img, mask, pad_r, pad_c)
% function [xx_fill, yy_fill, img_fill, mask_fill] = piv_fill(...
%   xx, yy, img, mask, pad_r, pad_c)
% 
% Fill non-sand regions of the input image by mirroring boundary pixels
%
% Arguments:
%   xx, yy: Vectors, coordinate axes for image columns, rows
%   img: 3D matrix, rectified/cropped RGB image
%   mask: 2D matrix, logical mask of sand pixels in img
%   pad_r, pad_c: Scalar integers, number of pixels of padding to add in
%       row and column directions
% 
% Outputs:
%   xx_fill, yy_fill: Vectors, coordinate axes for filled image columns, rows
%   img_fill: 3D matrix, filled and padded RGB image
%   mask_fill: 2D matrix, padded image mask matching (2D) size of img_fill
% %

% note: x-direction is padded with zeros and left masked. There is no
%   sane way to pad regions where sand grains enter or leave the frame.

% sanity checks
narginchk(6, 6);
validateattributes(xx, {'numeric'}, {'vector'});
validateattributes(yy, {'numeric'}, {'vector'});
validateattributes(img, {'numeric'}, {'3d', 'size', [numel(yy), numel(xx), 3]});
validateattributes(mask, {'logical'}, {'2d', 'size', [numel(yy), numel(xx)]});
validateattributes(pad_r, {'numeric'}, {'integer', 'nonnegative'});
validateattributes(pad_c, {'numeric'}, {'integer', 'nonnegative'});

fprintf('%s: fill and pad image and its mask by mirroring sand pixels\n', ...
    mfilename); 

% create padded arrays
img_fill = padarray(img, [pad_r, pad_c, 0], NaN, 'both');
mask_fill = padarray(mask, [pad_r, pad_c], 0, 'both');

% extend coordinate vectors
dy = yy(2) - yy(1);
yy_fill = [yy(1) - dy*(pad_r:-1:1), yy, yy(end) + dy*(1:pad_r)];

dx = xx(2) - xx(1);
xx_fill = [xx(1) - dx*(pad_c:-1:1), xx, xx(end) + dx*(1:pad_c)];

% temporarily convert the image to double
img_fill = double(img_fill);
    
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

% create (repeatable) random sawtooth pattern for mirror offsets
% note: using a constant seed to make this repeatable 
rng(204808, 'twister');  % set seed and generator type
max_depth = 20;

offset_idx = [];
while length(offset_idx) < size(img_fill, 1)
    depth = randi([0, max_depth]);
    offset_idx = [offset_idx, 1:depth, (depth-1):-1:2];  %#ok!
end

% FIXME: fails when layer is too thin?

% smooth the upper boundary line
% note: lower boundary is *not* smooth due to the presence of the 
%   metal support in the image, leave it alone
smooth_num_pts = 10;
top_bnd = smooth(top_row, smooth_num_pts/length(top_row), 'lowess')';
bot_bnd = bot_row;

% compute row index for fill region using offsets from the boundary
[nr, nc] = size(mask_fill);
[mirror_jj, mirror_ii] = meshgrid(1:nc, 1:nr);
for jj = 1:nc
    
    % fill above wedge
    if ~isnan(top_row(jj))
        fill_idx = top_row(jj):nr;  % portion of column to fill
        fill_values = mirror_ii(fill_idx, jj);
        fill_values = fill_values(offset_idx(1:length(fill_values)));  % apply sawtooth
        mirror_ii(fill_idx, jj) = top_bnd(jj) + top_bnd(jj) - fill_values;
    end
    
    % fill below wedge
    if ~isnan(bot_row(jj))
        fill_idx = bot_row(jj):-1:1;  % portion of column to fill
        fill_values = mirror_ii(fill_idx, jj);
        fill_values = fill_values(offset_idx(1:length(fill_values)));  % apply sawtooth
        mirror_ii(fill_idx, jj) = bot_bnd(jj) + bot_bnd(jj) - fill_values;
    end
end

% apply reflection by interpolating
% note: pixels within the mask are left unchanged, as desired
[original_jj, original_ii] = meshgrid(1:nc, 1:nr);
for cc = 1:3
    img_fill(:, :, cc) = interp2(...
        original_jj, original_ii, img_fill(:, :, cc), mirror_jj, mirror_ii, 'cubic');
end

% revert image to byte
img_fill = uint8(img_fill);

% sanity checks
assert(length(yy_fill) == size(img_fill, 1), 'Padded dimensions do not match padded image');
assert(length(xx_fill) == size(img_fill, 2), 'Padded dimensions do not match padded image');
assert(size(img_fill, 3) == 3, 'Padded image does not appear to be RGB');
