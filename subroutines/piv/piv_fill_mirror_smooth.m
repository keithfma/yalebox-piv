function [xx_fill, yy_fill, img_fill] = piv_fill_mirror_smooth(xx, yy, img, mask, pad_width)
% function [xx_fill, yy_fill, img_fill] = piv_fill_mirror_smooth(xx, yy, img, mask, pad_width)
% 
% Fill non-sand regions of the input image by mirroring boundary pixels,
% uses a weighted (gaussian) mean of boundary pixels to reduce sensitivity
% to the masking algorithm 
%
% Arguments:
%   xx, yy: Vectors, coordinate axes for image columns, rows
%   img: 3D matrix, rectified/cropped RGB image
%   mask: 2D matrix, logical mask of sand pixels in img
%   pad_width: Scalar integer, number of pixels of padding to add
% 
% Outputs:
%   xx_fill, yy_fill: Vectors, coordinate axes for filled image columns, rows
%   img_fill: 3D matrix, filled and padded RGB image
% %

% create padded arrays
img_fill = padarray(img, [pad_width, pad_width, 0], NaN, 'both');
mask_fill = padarray(mask, [pad_width, pad_width], 0, 'both'); % not returned

% create extended coordinate vectors
dx = xx(2) - xx(1);
xx_fill = [xx(1) - dx*(pad_width:-1:1), xx, xx(end) + dx*(1:pad_width)];

dy = yy(2) - yy(1);
yy_fill = [yy(1) - dy*(pad_width:-1:1), yy, yy(end) + dy*(1:pad_width)];

% create gaussian weight array
% note: kernel size must be odd
weight = fspecial('gaussian', 7, 1); 

% temporarily convert the image to double
img_fill = double(img_fill);

% fill upwards and downwards
% note: expand mask so that next step expands filled pixels to sides
for jj = 1:length(xx_fill)
    ii_bot = find(mask_fill(:, jj), 1, 'first');
    if ~isempty(ii_bot)
        for cc = 1:3
            img_fill(1:ii_bot, jj, cc) = get_value(img_fill, mask_fill, weight, ii_bot, jj, cc);
        end
        mask_fill(1:ii_bot, jj) = true;
    end
    ii_top = find(mask_fill(:, jj), 1, 'last');
    if ~isempty(ii_top)
        for cc = 1:3
            img_fill(ii_top:end, jj, cc) = get_value(img_fill, mask_fill, weight, ii_top, jj, cc);
        end
        mask_fill(ii_top:end, jj) = true;
    end
end

% note: filling the left and right boundaries will not work - the sand at
% these boundaries comes into / goes out of frame between images, and so
% there is no continuity. Leaving this step in for now to demonstrate this
% fact.

% fill leftwards and rightwards
% note: expand mask so that we can confirm all pixels were filled
for ii = 1:length(yy_fill)
    jj_left = find(mask_fill(ii, :), 1, 'first');
    if ~isempty(jj_left)
        for cc = 1:3
            img_fill(ii, 1:jj_left, cc) = img_fill(ii, jj_left, cc);
        end
        mask_fill(ii, 1:jj_left) = true;
    end
    jj_right = find(mask_fill(ii, :), 1, 'last');
    if jj_right > 4472
        keyboard
    end
    if ~isempty(jj_right)
        for cc = 1:3
            img_fill(ii, jj_right:end, cc) = img_fill(ii, jj_right, cc);
        end
        mask_fill(ii, jj_right:end) = true;
    end
end

% revert image to byte
img_fill = uint8(img_fill);

% sanity checks
assert(all(mask_fill(:)), 'Expected all pixels to be populated by fill routine');
assert(length(yy_fill) == size(img_fill, 1), 'Padded dimensions do not match padded image');
assert(length(xx_fill) == size(img_fill, 2), 'Padded dimensions do not match padded image');
assert(size(img_fill, 3) == 3, 'Padded image does not appear to be RGB');


function val = get_value(img, msk, wgt, ii, jj, cc)
% Return weighted image value at index (ii, jj, cc), taking care to handle
%   NaN-padded edges
% %
dd = (size(wgt,1)-1)/2;
msk_win = msk((ii-dd):(ii+dd), (jj-dd):(jj+dd)); 
wgt_win = wgt.*msk_win;
wgt_win = wgt_win./sum(wgt_win(:));
img_win = img((ii-dd):(ii+dd), (jj-dd):(jj+dd), cc);
val = sum(wgt_win(:).*img_win(:));

