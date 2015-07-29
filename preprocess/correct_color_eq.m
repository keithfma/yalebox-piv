function [valeq] = correct_color_eq(rgb, mask, width)

% convert image to 2D intensity matrix
hsv = rgb2hsv(rgb);
val = hsv(:,:,3);
[nrow, ncol] = size(val);

% equalize histogram, column-by-column
valeq = zeros([nrow, ncol]);
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
    valeq(mask_col) = interp1(old(2:end), new(2:end), val(mask_col));
end

% add a tiny offset so that no sand pixels are exactly zero
valeq(mask & valeq==0) = valeq(mask & valeq==0)+1e-5;