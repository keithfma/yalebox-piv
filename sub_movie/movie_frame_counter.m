function img_a = movie_frame_counter(img, x, y, count, pos, fsize, fclr, bclr, bopac)
% function img_a = movie_frame_counter(img, x, y, count, pos, fsize, fclr, bclr, bopac)
%
% Add an integer counter to the input image using the MATLAB insertText
% function.
%
% Arguments:
%
%   img_a = 3D matrix, the output annotated image.
%
%   img = 2D or 3D matrix, the input image to be annotated.
%
%   x, y = Vector, double,  coordinate vectors for the input image
%
%   count = Scalar integer, counter value to be added to the image.
%
%   pos = 2-element vector, double, location of the counter center in world
%       coordinates
%
%   fsize = Scalar, numeric, font size for title text.
%
%   fclr = MATLAB color definintion, can be a valid string or RGB triplet,
%       font color for title text.
%
%   bclr = MATLAB color definintion, can be a valid string or RGB triplet,
%       font color for text box.
%
%   bopac = Scalar double, normlized opacity of the text box, in the range
%       [0,1].
%
% Keith Ma, August 2015

% get dimensions
nx = numel(x);
ny = numel(y);

% convert position to pixel coords
pos(1) = interp1(x, 1:nx, pos(1), 'linear', 'extrap');
pos(2) = interp1(y, 1:ny, pos(2), 'linear', 'extrap');

% insert counter
img_a = insertText(img, pos, sprintf('step = %i', count), ...
    'FontSize', fsize, ...
    'TextColor', fclr, ...
    'BoxColor', bclr, ...
    'BoxOpacity', bopac, ...
    'AnchorPoint', 'Center'); 
