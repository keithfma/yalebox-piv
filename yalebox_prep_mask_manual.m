function mask = yalebox_prep_mask_manual(rgb)
% function mask = yalebox_prep_mask_manual(rgb)
%
% Interactively create a mask to black out region(s) of an image.
%
% Arguments:
%
%   rgb = 3D matrix, uint8, a 24-bit "Truecolor" RGB image, as read into
%       MATLAB with imread()
%
%   mask = 2D matrix, logical, true where there is sand and false
%       elsewhere.
%
% Keith Ma, July 2015

% check for sane arguments
assert(isa(rgb, 'uint8') && size(rgb,3) == 3, ...
    'img_rgb is not a 24-bit RGB image');

% initialize mask and figure
mask = false(size(rgb,1), size(rgb,2));
rgbm = rgb;
figure()
subplot(2,1,1)
imshow(rgb); title('original');
subplot(2,1,2); title('masked'); 

% add polygons to mask interactively
while 1
    
    % create new mask 
    subplot(2,1,1)
    p = impoly('Closed', true);
    this_mask = createMask(p);
    
    % apply mask and display
    rgbm(repmat(this_mask, 1, 1, 3)) = 0;
    subplot(2,1,2)
    imshow(rgbm);
    
    % query user
    choice = input('Keep (1) or discard (anything else): ');
    if choice == 1
        mask = mask | this_mask;
    end
    choice = input('Finished (1) or add another (anything else): ');
    if choice == 1
        break
    end
    
end
    
% invert mask
mask = ~mask;
