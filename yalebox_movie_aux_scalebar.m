function aimg = yalebox_movie_aux_scalebar(img, x, y, ulc, wid, ht, bclr, fclr, fsize)
%
% Add scalebar annotation to the input image using MATLABs insertShape and
% insertText functions.
%
% ulc = upper left corner, [x, y]

% get dimensions
nx = numel(x);
ny = numel(y);
dx = abs(x(1)-x(2));
dy = abs(y(1)-y(2));

% get position matrix for bar rectangle
ulx = interp1(x, 1:nx, ul(1), 'linear', 'extrap');
uly = interp1(y, 1:ny, ul(2), 'linear', 'extrap');



