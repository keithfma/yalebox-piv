function [h] = pcolorc(x,y,c)
% function [h] = pcolorc(x,y,c)
%
% Tweaked version of the MATLAB routine pcolor.  Creates a pseudocolor plot
% of the data c(x,y) with x,y in the center of each pixel and no rows or
% columns dropped.
%
% Keith Ma, September 2012

% get spacing and direction of coordinate grids
dx = x(1,2)-x(1,1);
xPositiveRight = dx>0; 
dx = abs(dx);

dy = y(2,1)-y(1,1);
yPositiveDown = dy>0; % this would be true for a standard meshgrid matrix
dy = abs(dy);

% flip matrices to a "standard" orientation
if ~xPositiveRight
    x = fliplr(x);
    c = fliplr(c);
end

if ~yPositiveDown
    y = flipud(y);
    c = flipud(c);
end

% add dummy row
x = [x; x(end,:)];
y = [y; (y(end,1)+dy)*ones(1,size(y,2))];
c = [c; nan(1,size(y,2))];

% add dummy column
x = [x, (x(1,end)+dx)*ones(size(x,1),1)];
y = [y, y(:,end)];
c = [c, nan(size(x,1),1)];

% shift coordinate grids to register properly
x = x-dx/2;
y = y-dy/2;

% plot using pcolor
h = pcolor(x,y,c);
