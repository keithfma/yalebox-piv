function [] = correct_color_eq(rgb, mask, width)

% convert image to 2D intensity matrix
hsv = rgb2hsv(rgb);
val = hsv(:,:,3);
[nrow, ncol] = size(val);
 

% equalize histogram, column-by-column
valeq = zeros([nrow, ncol]);
for i = 1:ncol
    
    ind = (i-width):(i+width);
    ind(ind<=0) = abs(ind(ind<=0))+2;
    ind(ind>ncol) = 2*ncol-ind(ind>ncol);
    
     val_win = val(:, ind);
     mask_win = mask(:, ind);
     eql = histeq(val_win(mask_win));
     val_win(mask_win) = eql;
     valeq(:,i) = val_win(:, width+1);
end
keyboard
    