function [xx_fill, yy_fill, img_fill, mask_fill] = piv_fill(...
    xx, yy, img, mask, pad_width, method)
% function [xx_fill, yy_fill, img_fill, mask_fill] = piv_fill(...
%   xx, yy, img, mask, pad_width, method)
% 
% Fill non-sand regions of the input image by extending boundary pixels
%
% Arguments:
%   xx, yy: Vectors, coordinate axes for image columns, rows
%   img: 3D matrix, rectified/cropped RGB image
%   mask: 2D matrix, logical mask of sand pixels in img
%   pad_width: Scalar integer, number of pixels of padding to add
%   method: (Optional) string, which boundary treatment option to apply,
%       default is 'extend_smooth'
% 
% Outputs:
%   xx_fill, yy_fill: Vectors, coordinate axes for filled image columns, rows
%   img_fill: 3D matrix, filled and padded RGB image
%   mask_fill: 2D matrix, padded image mask matching (2D) size of img_fill
% %

% FIXME: remove the method option once the dust has settled on image padding,
%   there should be one and only one fill method implemented here in v1.0

% set defaults
narginchk(5, 6);
if nargin == 5; method = 'extend_smooth'; end

fprintf('%s: apply boundary treatment "%s" to input image\n', ...
    mfilename, method); 

if strcmp(method, 'mask')
    % legacy, don't pad just apply mask to image 
    img_fill = img;
    img_fill(repmat(~mask, 1, 1, size(img_ti, 3))) = 0;
    xx_fill = xx;
    yy_fill = yy;
    mask_fill = mask;

elseif strcmp(method, 'mirror')
    % mirror boundary pixels
    [xx_fill, yy_fill, img_fill, mask_fill] = piv_fill_mirror(...
        xx, yy, img, mask, pad_width);
    
else
    error('bad value for argument "method": %s', method);
    
end


function [xx_fill, yy_fill, img_fill, mask_fill] = piv_fill_mirror(...
    xx, yy, img, mask, pad_width)
% function [xx_fill, yy_fill, img_fill, mask_fill] = piv_fill_mirror(...
%   xx, yy, img, mask, pad_width)
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
%   mask_fill: 2D matrix, padded image mask matching (2D) size of img_fill
% %

% create padded arrays
img_fill = padarray(img, [pad_width, 0, 0], NaN, 'both');
mask_fill = padarray(mask, [pad_width, 0], 0, 'both');

% extend coordinate vectors
dy = yy(2) - yy(1);
yy_fill = [yy(1) - dy*(pad_width:-1:1), yy, yy(end) + dy*(1:pad_width)];

xx_fill = xx; 

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

% build index array by reflecting until all indices are within sand
[nr, nc] = size(mask_fill);
bot_rows = repmat(bot_row, nr, 1);
top_rows = repmat(top_row, nr, 1);
[cols, rows] = meshgrid(1:nc, 1:nr);

while any(rows(:) > top_rows(:) | rows(:) < bot_rows(:))
    
    % reflect at top boundary
    above_top = rows > top_rows;
    rows(above_top) = 2*top_rows(above_top) - rows(above_top);  % same as t-(r-t)
    
    % reflect at bottom boundary
    below_bot = rows < bot_rows;
    rows(below_bot) = 2*bot_rows(below_bot) - rows(below_bot);  % same as b+(b-r)

end
idx = sub2ind([nr, nc], rows, cols);

% apply reflection
for cc = 1:3
    band = img_fill(:, :, cc);
    band(:) = band(idx);
    img_fill(:, :, cc) = band;
end

% revert image to byte
img_fill = uint8(img_fill);

% sanity checks
assert(length(yy_fill) == size(img_fill, 1), 'Padded dimensions do not match padded image');
assert(length(xx_fill) == size(img_fill, 2), 'Padded dimensions do not match padded image');
assert(size(img_fill, 3) == 3, 'Padded image does not appear to be RGB');
