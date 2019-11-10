function [img_fill, mask_fill] = prep_fill_mirror(img, mask)
% function img_fill = prep_fill(img, mask)
% 
% Fill non-sand regions of the input image by mirroring across sand boundary
%
% Arguments:
%   img: 2D matrix, rectified/cropped, and padded grayscale image
%   mask: 2D matrix, logical mask of sand pixels in img
%   method: String, name of method to use for filling pixels, current options are: "none", default
%       is "none"
%   show: Logical, set true to display results for visual debugging, default is false
% 
% Outputs:
%   img_fill: 2D matrix, filled grayscale image
%   mask_fill: 2D matrix, image mask matching size of img_fill
% %

% note: x-direction is padded with zeros and left masked. There is no
%   sane way to pad regions where sand grains enter or leave the frame.

% find row index of the top and bottom boundaries
[nr, nc] = size(img);
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
img_fill = img;
[jj, ii] = meshgrid(1:nc, 1:nr); 
img_fill(:, :) = interp2(jj, ii, img_fill(:, :), cols, rows, 'bilinear');  % FIXME: higher order interpolants introduce NaNs

% black out regions with no sand at all
mask_fill = ~has_sand;
img_fill(~has_sand) = NaN;



