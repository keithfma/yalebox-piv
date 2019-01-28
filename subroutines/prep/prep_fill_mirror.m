function [xx_fill, yy_fill, img_fill] = prep_fill_mirror(xx, yy, img, mask, pad_width)
% function img = prep_fill_mirror(img, mask, pad_width)
% 
% Fill non-sand regions of the input image by mirroring boundary pixels
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
img_fill = padarray(img, [pad_width, pad_width, 0], 0, 'both');
mask_fill = padarray(mask, [pad_width, pad_width], 0, 'both'); % not returned

% create extended coordinate vectors
dx = xx(2) - xx(1);
xx_fill = [xx(1) - dx*(pad_width:-1:1), xx, xx(end) + dx*(1:pad_width)];

dy = yy(2) - yy(1);
yy_fill = [yy(1) - dy*(pad_width:-1:1), yy, yy(end) + dy*(1:pad_width)];

% fill upwards and downwards
% note: expand mask so that next step expands filled pixels to sides
for jj = 1:length(xx_fill)
    ii_bot = find(mask_fill(:, jj), 1, 'first');
    if ~isempty(ii_bot)
        for cc = 1:3
            img_fill(1:ii_bot, jj, cc) = img_fill(ii_bot, jj, cc);
        end
        mask_fill(1:ii_bot, jj) = true;
    end
    ii_top = find(mask_fill(:, jj), 1, 'last');
    if ~isempty(ii_top)
        for cc = 1:3
            img_fill(ii_top:end, jj, cc) = img_fill(ii_top, jj, cc);
        end
        mask_fill(ii_top:end, jj) = true;
    end
end

% TODO: this is clearly the wrong approach at LR boundaries...
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

% sanity checks
assert(all(mask_fill(:)), 'Expected all pixels to be populated by fill routine');
assert(length(yy_fill) == size(img_fill, 1), 'Padded dimensions do not match padded image');
assert(length(xx_fill) == size(img_fill, 2), 'Padded dimensions do not match padded image');
assert(size(img_fill, 3) == 3, 'Padded image does not appear to be RGB');


