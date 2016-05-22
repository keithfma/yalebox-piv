function img_a = movie_frame_title(img, str, fsize, fclr, bclr, bopac)
% function img_a = movie_frame_title(img, str, fsize, fclr, bclr, bopac)
%
% Add title annotation to the top side of the input image using the MATLAB
% insertText function. Assumes that the top of the image is at row == 1.
%
% Arguments:
%
%   img_a = 3D matrix, the output annotated image.
%
%   img = 2D or 3D matrix, the input image to be annotated.
%
%   str = String, title text to be added to the image.
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

% get postion
pos = [size(img,2)/2, 1];

% annotate image
img_a = insertText(img, pos, str, ...
    'FontSize', fsize, ...
    'TextColor', fclr, ...
    'BoxColor', bclr, ...
    'BoxOpacity', bopac, ...
    'AnchorPoint', 'CenterTop');      
