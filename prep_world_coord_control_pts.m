function woco = prep_world_coord_control_pts(woco_image, show)
% function woco = prep_world_coord_control_pts(woco_image, show)
%
% Interactively extract control points defining the world coordinate image.
%
% Arguments:
%
%   woco_image = String, filename of the world coordinate grid image,
%       required.
%
%   show = (Optional) Scalar, logical flag, default = true
%
%   woco.xw, woco.yw = 1D vectors, control point world coordinate x- and
%       y-position in meters
%
%   woco.xp, woco.yp = 1D vectors, control point image coordinate x- and
%       y-position in pixels
%
% Keith Ma

% set defaults
if nargin == 1; 
    show = true; 
end

% check for sane arguments
narginchk(1, 2);
validateattributes(woco_image, {'char'}, {'vector'});
validateattributes(show, {'numeric', 'logical'}, {'scalar'});

% display coordinate grid image
im = imread(woco_image);
imshow(im)
axis on;
hold on
ax = gca;

% select npts control points interactively
pts = impoint(ax);
while 1
    switch questdlg('Add another control point?', 'prep_get_world_coord',  'Yes', 'No', 'No'); 
        case 'Yes'
            % do nothing
        case 'No'
            disp('n')
            break
    end
    pts(end+1) = impoint(ax); %#ok
end

% enter control point positions in world coordinates
npts = length(pts);
xw = zeros(npts, 1);
yw = zeros(npts, 1);
ii = 1;
while 1
    %...show current point
    tmp = pts(ii).getPosition();
    hpt = plot(ax, tmp(1), tmp(2), 'Color', 'r', 'MarkerSize', 20, 'Marker', 'o', 'LineStyle', 'None');
    xlabel(sprintf('current point x = %.2f, y = %.2f', xw(ii), yw(ii)));
    %...enter world coordinates
    try
        tmp = inputdlg({'X', 'Y'}, 'Input world coordinates in meters as X, Y', 1, ...
            {sprintf('%.2f', xw(ii)), sprintf('%.2f', yw(ii))} );
        xw(ii) = str2double(tmp{1});
        yw(ii) = str2double(tmp{2});
    catch
        %...repeat
        fprintf('invalid response, try again');
        delete(hpt);
        continue
    end
    %...done?
    if ii == npts
        switch questdlg('Satisfied with world coordinates for all control points?', 'prep_get_world_coord',  'Yes', 'No', 'No');
            case 'Yes'
                delete(hpt);
                break
            case 'No'
                % do nothing
        end
    end
    %...next
    ii = 1+mod(ii, npts);
    delete(hpt);
end


% opportunity to fine tune point location
while 1
    if input('Enter (1) when all control point locations are correct: ') == 1
        break
    end
end

% record control point positions in image coordinates (col, row)
xp = zeros(npts, 1);
yp = zeros(npts, 1);
for ii = 1:npts
    tmp = pts(ii).getPosition();
    xp(ii) = tmp(1);
    yp(ii) = tmp(2);
end
 
% prepare output struct
woco = struct('xw', xw, 'yw', yw, 'xp', xp, 'yp', yp);

% (optional) show image with annotated control points
if show
    figure()
    imshow(im);
    axis on
    title('Control points')
    xlabel('X [pixels]');
    ylabel('Y [pixels]');
    hold on
    for ii = 1:npts
        plot(woco.xp(ii), woco.yp(ii), ...
            'Marker', 'x', 'MarkerSize', 10, 'Color', 'r', 'LineStyle', 'None');
        text(woco.xp(ii), woco.yp(ii), sprintf('(%.2f, %.2f)', woco.xw(ii), woco.yw(ii)), ...
            'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end