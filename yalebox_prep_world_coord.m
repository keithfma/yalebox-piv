function [x, y] = yalebox_prep_world_coord(woco_image)
%
%
% Arguments:
%
%   woco_image = String, filename of the world coordinate grid image.
%
%   x = 1D vector, length == number of columns in coordinage grid image,
%       double, world coordinate x-position in meters 
%
%   y = 1D vector, length == number of rows in coordinage grid image,
%       double, world coordinate y-position in meters 
%
% Keith Ma, July 2015

% check for sane arguments

% display coordinate grid image
im = imread(woco_image);
imshow(im)
hold on

% select control points interactively - minimum of 4 points, no maximum
ctrl_pt_obj = [impoint(), impoint(), impoint()];
while 1
    ctrl_pt_obj(end+1) = impoint();
    if input('Enter (1) to stop adding points: ') == 1;
        break
    end
end
npts = numel(ctrl_pt_obj);

% get control point positions in image coordinates (col, row)
ctrl_pt_image = nan(npts, 2);
for i = 1:npts
    ctrl_pt_image(i,:) = ctrl_pt_obj(i).getPosition();
    delete(ctrl_pt_obj(i));
end

% review control points and enter positions
ctrl_pt_world = nan(npts, 2);
i = 1;
while i <= npts
    
    % request x, y position from user
    this = plot(ctrl_pt_image(i,1), ctrl_pt_image(i,2), '*r');
    try 
        tmp = inputdlg({'X', 'Y'}, 'Input world coordinates in meters as X, Y');
        delete(this);
        ctrl_pt_world(i, :) = str2num(cell2mat(tmp));
    catch
        continue
    end
    
    % next point
    plot(ctrl_pt_image(i,1), ctrl_pt_image(i,2), '*k');
    text(ctrl_pt_image(i,1), ctrl_pt_image(i,2), ...
        sprintf('%.2f, %.2f', ctrl_pt_world(i,1), ctrl_pt_world(i,2)),...
        'VerticalAlignment', 'Bottom',...
        'HorizontalAlignment', 'Center');
    i = i+1;
end


% best-fit scale and offset parameters


% create coordinate vectors


keyboard