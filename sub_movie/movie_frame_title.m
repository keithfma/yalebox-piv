function img_a = movie_frame_title(img, str, fsize, fclr, bclr, bopac, loc)
% function img_a = movie_frame_title(img, str, fsize, fclr, bclr, bopac, loc)
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
%   loc = String, horizontal position of title in window, i.e. 'left', 'center', 'right'
%
% % Keith Ma

% get postion and anchor
if strcmp(loc, 'left')    
    pos = [1, 1];
    anchor = 'LeftTop';
elseif strcmp(loc, 'center')
    pos = [size(img,2)/2, 1];
    anchor = 'CenterTop';
elseif strcmp(loc, 'right')
    pos = [size(img,2), 1];
    anchor = 'RightTop';
else
    error('Bad value for argument "loc"');
end
    


% annotate image
img_a = insertText(img, pos, str, ...
    'FontSize', fsize, ...
    'TextColor', fclr, ...
    'BoxColor', bclr, ...
    'BoxOpacity', bopac, ...
    'AnchorPoint', anchor);      
