function [out, prm] = yalebox_prep_intensity(rgb, mask, width)
% function [out, prm] = yalebox_prep_intensity(rgb, mask, width)
%
% Compute and apply local histogram-equalization color corrections.
% Corrections are computed for each column based on the histogram for sand
% pixels within "width" columns. This sliding-window approach corrects for
% lighting gradients in the x-direction (but not the y-direction). For
% details on the equalization algorithm, see:
%
% http://fourier.eng.hmc.edu/e161/lectures/contrast_transform/node2.html
%
% All arguments are required - default values are used where inputs are []
%
% Arguments:
%
%   rgb = 3D matrix, uint8, a 24-bit "Truecolor" RGB image, as read into
%       MATLAB with imread()
%
%   mask = 2D matrix, double, TRUE where there is sand, FALSE elsewhere
%
%   width = Scalar, double, half-width of the sliding window for quantile
%       calculations, in (whole) pixels, default = 50; 
%
%   out = 2D matrix, size(mask), double, normalized intensity image derived from
%     rgb, histogram is approximately uniform
%
%   prm = Struct,  contains all of the parameters (input and internally
%       computed), included so that default values can be recovered.
%       Structure member variables are: width, tiny
%
% Keith Ma, July 2015

% set defaults
if isempty(width); width = 50; end

% check for sane inputs
assert(isa(rgb, 'uint8') && size(rgb,3) == 3, ...
     'rgb is not a Truecolor (24-bit) RGB image');
assert(isa(mask, 'logical'), ...
    'mask is not of type logical ');
assert(size(mask,1) == size(rgb, 1) && size(mask,2) == size(rgb, 2), ...
    'mask and rgb dimensions do not match');
assert(numel(width) == 1 && mod(width,1) == 0, ...
    'width is not a scalar integer');

% convert image to 2D intensity matrix
hsv = rgb2hsv(rgb);
val = hsv(:,:,3);
[nrow, ncol] = size(val);

% equalize histogram, column-by-column
out = zeros([nrow, ncol]);
for i = 1:ncol
    
    % window indices, with symetric padding
    ind = (i-width):(i+width);
    ind(ind<=0) = abs(ind(ind<=0))+2;
    ind(ind>ncol) = 2*ncol-ind(ind>ncol);
    
    % mask for sand pixels in col and in window
    mask_col = false(size(mask));
    mask_col(:,i) = mask(:,i);
    mask_win = false(size(mask));
    mask_win(:,ind) = mask(:,ind);
    
    % lookup table (empirical cdf for sand pixels in window)
    [new, old] = ecdf(val(mask_win));
    
    % correct using lookup table
    out(mask_col) = interp1(old(2:end), new(2:end), val(mask_col));
end

% add a tiny offset so that no sand pixels are exactly zero
tiny = 1e-5;
out(mask & out==0) = out(mask & out==0)+tiny;

% prepare parameter struct
prm.width = width;
prm.tiny = tiny;
