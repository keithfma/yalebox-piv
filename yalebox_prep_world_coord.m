function [x, y] = yalebox_prep_world_coord(woco_image, npts)
%
%
% Arguments:
%
%   woco_image = String, filename of the world coordinate grid image.
%
%   npts = Integer, number of control points to define 
%
%   x = 1D vector, length == number of columns in coordinage grid image,
%       double, world coordinate x-position in meters 
%
%   y = 1D vector, length == number of rows in coordinage grid image,
%       double, world coordinate y-position in meters 
%
% Keith Ma, July 2015

% check for sane arguments
% npts >= 4

% display coordinate grid image
im = imread(woco_image);
imshow(im)
hold on

% select npts control points interactively
ctrl_pt_obj = impoint();
for i = 2:npts
    ctrl_pt_obj(end+1) = impoint(); %#ok
end
while 1
    if input('Enter (1) when all control points are correct: ') == 1
        break
    end
end

% get control point positions in image coordinates (col, row)
ctrl_pt_image = nan(npts, 2);
for i = 1:npts
    ctrl_pt_image(i,:) = ctrl_pt_obj(i).getPosition();
    delete(ctrl_pt_obj(i));
end

% review control points and enter positions in world coordinate (x, y)
ctrl_pt_world = nan(npts, 2);
i = 1;
while i <= npts
    
    % request x, y position from user
    this = plot(ctrl_pt_image(i,1), ctrl_pt_image(i,2), '*r');
    try 
        tmp = inputdlg({'X', 'Y'}, 'Input world coordinates in meters as X, Y');
        delete(this);
        ctrl_pt_world(i, 1) = str2num(tmp{1});
        ctrl_pt_world(i, 2) = str2num(tmp{2});
    catch err
        if isempty(tmp) 
            return
        else
            continue
        end
    end
    
    % next point
    plot(ctrl_pt_image(i,1), ctrl_pt_image(i,2), '*k');
    text(ctrl_pt_image(i,1), ctrl_pt_image(i,2), ...
        sprintf('%.2f, %.2f', ctrl_pt_world(i,1), ctrl_pt_world(i,2)),...
        'VerticalAlignment', 'Bottom',...
        'HorizontalAlignment', 'Center');
    i = i+1;
end

% best-fit for transformation equation world = m*pixel+b

% x-dir
param = [ctrl_pt_image(:,1), ones(npts,1)]\ctrl_pt_world(i,:);

keyboard
% % x-dir
% A = [ones(npts,1), xp(:), yp(:)];
% param = A\xw;
% xw0 = param(1);
% Jp2w(1,1) = param(2);
% Jp2w(1,2) = param(3);
% 
% % y-dir
% param = A\yw;
% yw0fit = param(1);
% Jp2w(2,1) = param(2);
% Jp2w(2,2) = param(3);


% create coordinate vectors

keyboard