function mask_manual = yalebox_prep_mask_manual(img)
% function mask_manual = yalebox_prep_mask_manual(img)
%
% Interactively create a mask to black out region(s) of an image.
%
% Arguments:
%
%   img = Matrix, image data that can be displayed by imshow
%
%   mask_manual = 2D matrix, logical, true where there is sand and false
%       elsewhere.
%
% Keith Ma, July 2015

% check for sane arguments
validateattributes(img, {'numeric'}, {'3d', '2d'}, ...
    'yalebox_prep_mask_manual', 'img');

% init image and mask
imgm = img;
mask = false(size(img,1), size(img,2));

% init figure
figure()
subplot(2,1,1)
imshow(img); title('original');
subplot(2,1,2); title('masked'); 

% add polygons to mask interactively
while 1
    
    % create new mask 
    subplot(2,1,1)
    p = impoly('Closed', true);
    this_mask = createMask(p);
    
    % apply mask and display
    imgm(repmat(this_mask, 1, 1, size(img, 3))) = 0;
    subplot(2,1,2)
    imshow(imgm);
    
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
    
% invert mask and rename
mask_manual = ~mask;
