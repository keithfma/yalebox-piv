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
im = imread(woco_image);

% select control points interactively - minimum of 4 points, no maximum
imshow(im)
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

% MAKE THIS A FINITE WHILE LOOP, ONLY INCREMENT IF THERE IS NO ERROR

% review control points and enter positions
ctrl_pt_world = nan(npts, 2);
hold on
for i = 1:npts
    
    % request x, y position from user
    this = plot(ctrl_pt_image(i,1), ctrl_pt_image(i,2), '*r');
    tmp = inputdlg({'X', 'Y'}, 'Input world coordinates in meters as X, Y');
    ctrl_pt_world(i, :) = str2num(cell2mat(tmp));
    
    % update display
    delete(this);
    plot(ctrl_pt_image(i,1), ctrl_pt_image(i,2), '*k');
    text(ctrl_pt_image(i,1), ctrl_pt_image(i,2), ...
        sprintf('%.2f, %.2f', ctrl_pt_world(i,1), ctrl_pt_world(i,2)),...
        'VerticalAlignment', 'Bottom',...
        'HorizontalAlignment', 'Center');
end


% best-fit scale and offset parameters


% create coordinate vectors


keyboard