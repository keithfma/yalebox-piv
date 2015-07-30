function [x, y, prm] = yalebox_prep_world_coord(woco_image, npts, show)
% function [x, y, prm] = yalebox_prep_world_coord(woco_image, npts, show)
%
% Compute the best-fit cartesian coordinate system for the coordinate grid
% image "woco_image". As a result of previous preprocessing steps (i.e.
% distortion correction and cropping), the coordinate system is neither
% rotated or skew. This means that the coordinate transformation can be
% descrbed with four free parameters, and the x- and y-direction transforms
% are independent. Specifically, the transformation is:
%
%   x_world = x_pixel*x_scale + x_offset
%   y_world = y_pixel*y_scale + y_offset
%
% Arguments:
%
%   woco_image = String, filename of the world coordinate grid image,
%       required.
%
%   npts = Integer, number of control points to define, must be >= 4 to
%       provide sufficient infomation to fit four coordinate transformation
%       parameters, default = 4
%
%   show = True to display image in world coordinates, default = true
%
%   x = 1D vector, length == number of columns in coordinage grid image,
%       double, world coordinate x-position in meters 
%
%   y = 1D vector, length == number of rows in coordinage grid image,
%       double, world coordinate y-position in meters 
%
%   prm = Struct, contains all of the parameters (input and internally
%       computed), included so that default values can be recovered.
%       Structure member variables are: woco_image, npts, x_scale,
%       x_offset, y_scale, y_offset
%
% Keith Ma, July 2015

% set defaults
if isempty(npts); npts = 4; end
if isempty(show); show = true; end

% check for sane arguments
assert(isa(woco_image, 'char') && exist(woco_image, 'file') == 2, ...
    'woco_image is not a valid filename');
assert(numel(npts) == 1 && mod(npts, 1) == 0 && npts >= 4, ...
    'npts is not a scalar integer greater than 4.');
assert(numel(show) == 1 && (show == 1 || show == 0), ...
    'show is not a scalar logical flag.');

% display coordinate grid image
im = imread(woco_image);
imshow(im)
hold on

% select npts control points interactively
pts = impoint();
for i = 2:npts
    pts(end+1) = impoint(); %#ok
end
while 1
    if input('Enter (1) when all control points are correct: ') == 1
        break
    end
end

% get control point positions in image coordinates (col, row)
xp = nan(npts,1);
yp = nan(npts,1);
for i = 1:npts
    tmp = pts(i).getPosition();
    xp(i) = tmp(1);
    yp(i) = tmp(2);
    delete(pts(i));
end

% review control points and enter positions in world coordinate (x, y)
xw = nan(npts,1);
yw = nan(npts,1);
i = 1;
while i <= npts
    
    % request x, y position from user
    this = plot(xp(i), yp(i), '*r');
    try 
        tmp = inputdlg({'X', 'Y'}, 'Input world coordinates in meters as X, Y');
        delete(this);
        xw(i) = str2num(tmp{1});
        yw(i) = str2num(tmp{2});
    catch
        if isempty(tmp) 
            return
        else
            continue
        end
    end
    
    % next point
    plot(xp(i), yp(i), '*k');
    text(xp(i), yp(i), ...
        sprintf('%.2f, %.2f', xw(i), yw(i)),...
        'VerticalAlignment', 'Bottom',...
        'HorizontalAlignment', 'Center');
    i = i+1;
end

% best-fit for transformation equation world = scale*pixel+offset
tmp = [xp, ones(npts,1)]\xw;
x_scale = tmp(1);
x_offset = tmp(2);

tmp = [yp, ones(npts,1)]\yw;
y_scale = tmp(1);
y_offset = tmp(2);

% create coordinate vectors
x = (1:size(im,2))*x_scale+x_offset;
y = (1:size(im,1))*y_scale+y_offset;

% (optional) show image with best-fit world coordinates
if show
    figure()
    imagesc(x,y,rgb2gray(im));
    colormap(gray);
    grid on
    title('Best-fit world coordinates')
    xlabel('X, meters');
    ylabel('Y, meters');
end

% prepare parameter struct
prm.woco_image = woco_image;
prm.npts = npts;
prm.x_scale = x_scale;
prm.x_offset = x_offset;
prm.y_scale = y_scale;
prm.y_offset = y_offset;

