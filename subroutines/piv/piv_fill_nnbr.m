function [xx_fill, yy_fill, img_fill] = piv_fill_nnbr(xx, yy, img, mask, pad_width)
% function [xx_fill, yy_fill, img_fill] = piv_fill_nnbr(xx, yy, img, mask, pad_width)
% 
% Fill non-sand regions of the input image by nearest-neighbor
% extrapolation
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

% find nearest neighbors outside sand region
[~, nbr_idx] = bwdist(mask_fill);
nbr_idx = nbr_idx(~mask_fill);

% fill non-sand region with nearest neighbor values
for cc = 1:3
    band = img_fill(:, :, cc);
    band(~mask_fill) = band(nbr_idx);
    img_fill(:, :, cc) = band;
end

% sanity checks
assert(length(yy_fill) == size(img_fill, 1), 'Padded dimensions do not match padded image');
assert(length(xx_fill) == size(img_fill, 2), 'Padded dimensions do not match padded image');
assert(size(img_fill, 3) == 3, 'Padded image does not appear to be RGB');

