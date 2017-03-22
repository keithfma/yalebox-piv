function img_a = movie_frame_scalebar(img, x, y, str, pos, clr, opac, fsize, fclr, bclr, bopac)
% function img_a = movie_frame_scalebar(img, x, y, str, pos, clr, opac, fsize, fclr, bclr, bopac)
%
% Add scalebar annotation to the input image using MATLABs insertShape and
% insertText functions.
%
% Arguments:
%
%   img_a = 3D matrix, the output annotated image.
%
%   img = 2D or 3D matrix, the input image to be annotated.
%
%   str = String, label fto be added above the scalebar.
%
%   pos = 4-element vector, scalebar position and dimensions is world
%       coordinates as [x_upper_left, y_upper_left, width, height]
%
%   clr = MATLAB color definintion, can be a valid string or RGB triplet,
%       scalebar color
%
%   opac = Scalar double, normalized opacity of the scalebar, in the range
%       [0,1].
%
%   fsize = Scalar, numeric, font size for label text.
%
%   fclr = MATLAB color definintion, can be a valid string or RGB triplet,
%       font color for label text.
%
%   bclr = MATLAB color definintion, can be a valid string or RGB triplet,
%       font color for containing box.
%
%   bopac = Scalar double, normlized opacity of the containing box, in the range
%       [0,1].
%
% Keith Ma, August 2015

% get dimensions
nx = numel(x);
ny = numel(y);
dx = abs(x(1)-x(2));
dy = abs(y(1)-y(2));

% convert position to pixel coords
pos(1) = interp1(x, 1:nx, pos(1), 'linear', 'extrap');
pos(2) = interp1(y, 1:ny, pos(2), 'linear', 'extrap');
pos(3) = pos(3)/dx;
pos(4) = pos(4)/dy;

% add text
img_a = insertText(img, [pos(1)+pos(3)/2, pos(2)+pos(4)/2], str, ...
    'FontSize', fsize, ...
    'TextColor', fclr, ...
    'BoxColor', bclr, ...
    'BoxOpacity', bopac, ...
    'AnchorPoint', 'CenterBottom');

% add scalebar 
img_a = insertShape(img_a, 'FilledRectangle', pos, ...
    'Color', clr, ...
    'Opacity', opac);