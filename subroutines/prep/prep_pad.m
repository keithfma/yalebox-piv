function [x_pad, y_pad, img_pad, mask_pad] = prep_pad(x_vec, y_vec, img, mask, pad_r, pad_c)
% function [x_pad, y_pad, img_pad, mask_pad] = prep_pad(x_vec, y_vec, img, mask, pad_r, pad_c)
% 
% Add pad to image, its mask, and its coordinate vectors
%
% Arguments:
%   x_vec, y_vec: Vectors, coordinate axes for image columns, rows
%   img: 2D matrix, rectified/cropped grayscale image
%   mask: 2D matrix, logical mask of sand pixels in img
%   pad_r, pad_c: Scalar integers, number of pixels of padding to add in
%       row and column directions
% 
% Outputs:
%   x_pad, y_pad: Vectors, coordinate axes for filled image columns, rows
%   img_pad: 2D matrix, filled and padded grayscale image
%   mask_pad: 2D matrix, padded image mask matching (2D) size of img_pad
% %

% note: x-direction is padded with zeros and left masked. There is no
%   sane way to pad regions where sand grains enter or leave the frame.

% sanity checks
narginchk(6, 6);
validateattributes(x_vec, {'numeric'}, {'vector'});
validateattributes(y_vec, {'numeric'}, {'vector'});
validateattributes(img, {'double', 'single'}, {'2d', 'size', [numel(y_vec), numel(x_vec)]});
validateattributes(mask, {'logical'}, {'2d', 'size', [numel(y_vec), numel(x_vec)]});
validateattributes(pad_r, {'numeric'}, {'integer', 'nonnegative'});
validateattributes(pad_c, {'numeric'}, {'integer', 'nonnegative'});

fprintf('%s: pad image, its mask, and its coordinate arrays by %d rows and %d columns\n', ...
    mfilename, pad_r, pad_c); 

% create padded arrays
img_pad = padarray(img, [pad_r, pad_c], NaN, 'both');  % fill is arbitrary, we use mask to find sand 
mask_pad = padarray(mask, [pad_r, pad_c], false, 'both');  

% extend coordinate vectors
dy = y_vec(2) - y_vec(1);
y_pad = [y_vec(1) - dy*(pad_r:-1:1), y_vec, y_vec(end) + dy*(1:pad_r)];

dx = x_vec(2) - x_vec(1);
x_pad = [x_vec(1) - dx*(pad_c:-1:1), x_vec, x_vec(end) + dx*(1:pad_c)];