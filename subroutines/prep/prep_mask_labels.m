function [labels, poly_out] = prep_mask_labels(rgb, poly_in, interactive)
% function [labels, polygons] = prep_mask_labels(rgb, polygons, interactive)
% 
% Interactive tool for defining pixel classes using ROI polygons. Select
% the region(s) that represent a single class (e.g. sand or other)
% 
% Arguments:
%   rgb: 3D, 24-bit color image
%
%  poly_in = initial ROI polygons for labeled image area, polygon vertices
%    as 2d array with points in columns and NaNs separating each polygon
%
%  interactive: logical scalar, set true to define/update label polygons
%    interactively, default is false
% 
% Returns: 
%  labels: integer array, linear indices into one layer of RGB for all
%    pixels inside the label polygons
% 
%  poly_out = final ROI polygons for labeled image area, polygon vertices
%    as 2d array with points in columns and NaNs separating each polygon
% % 

% set defaults
if nargin < 3; interactive = false; end

% TODO: sanity checks

if interactive
    % interactive polygon definition
    unpacked_poly_out = define_polygons(rgb, unpack_polygons(poly_in));
else
    % use input polygons without change
    unpacked_poly_out = unpack_polygons(poly_in);
end

mask = false(size(rgb, 1), size(rgb, 2));
for ii = 1:length(unpacked_poly_out)
    x_pts = unpacked_poly_out{ii}(:, 1);
    y_pts = unpacked_poly_out{ii}(:, 2);
    mask = mask | poly2mask(x_pts, y_pts, size(mask, 1), size(mask, 2));
end 
labels = find(mask);

poly_out = pack_polygons(unpacked_poly_out);


function outpoly = define_polygons(rgb, inpoly)
% Interactive GUI for defining/updating label ROI polygons
% 
% Arguments:
%   rgb: 3D, 24-bit color image
% 
%   inpoly = unpacked polygons, cell array or 2D matrices, each cell
%       contains vertices for one polygon with x-coord in row 1 and y-coord
%       in row 2
% 
% Returns:
%   outpoly = unpacked polygons, cell array or 2D matrices, each cell
%       contains vertices for one polygon with x-coord in row 1 and y-coord
%       in row 2
% %

% constants
margin = 0.025; % norm
buffer = 0.01; % norm
font = 14; % pt
instruct = {'+ Create a polygon: click "Add" button and then click to define vertices', ...
            '+ Edit a polygon: click andd vertices', ...
            '+ Delete existing: right click and select "Delete Polygon"', ...
            '+ Finish: click "Done"'};

img_left = margin;
img_bot = 0.1;
img_width = 1 - 2*margin;
img_height = 1 - (img_bot + 2*margin);
img_pos = [img_left, img_bot, img_width, img_height];

add_left = margin;
add_bot = margin;
add_width = 0.1;
add_height = 0.05;
add_pos = [add_left, add_bot, add_width, add_height];

done_left = add_left + add_width + buffer;
done_bot = add_bot;
done_width = add_width;
done_height = add_height;
done_pos = [done_left, done_bot, done_width, done_height];

help_left = done_left + done_width + buffer;
help_bot = done_bot;
help_width = done_width;
help_height = done_height;
help_pos = [help_left, help_bot, help_width, help_height];

% create UI
hf = figure('Units', 'Normalized', 'Outerposition', [0 0 1 1]);

hax = axes('Units', 'Normalized', 'Position', img_pos);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', add_pos, 'String', 'Add', 'FontSize', font, ...
    'Callback', @do_add);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', done_pos, 'String', 'Done', 'FontSize', font, ...
    'Callback', @(~,~)uiresume);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', help_pos, 'String', 'Help', 'FontSize', font, ...
    'Callback', @(~,~)msgbox(instruct));

% display image and input polygons 
imshow(rgb);
for ii = 1:length(inpoly)
    images.roi.Polygon(gca, 'Tag', 'ROI', 'Position', inpoly{ii});
end
hax.YDir = 'normal';

% wait here until user presses "Done"
uiwait(); 

% retreive and reformat output polygons
objs = findobj(hf, 'Tag', 'ROI');
outpoly = cell(size(objs));
for jj = 1:length(objs)
    outpoly{jj} = objs(jj).Position;
end

close(hf)


function [] = do_add(~, ~)
% Callback for add mask layer button
% %
drawpolygon('Tag', 'ROI');
uiwait();


function unpacked = unpack_polygons(packed)
% function unpacked = unpack_polygons(packed)
% 
% Arguments:
%   packed = polygon vertices as 2d array with points in columns and NaNs
%       separating each polygon
%
% Returns:
%   unpacked = cell array or 2D matrices, each cell contains vertices for
%       one polygon with x-coord in row 1 and y-coord in row 2
% % 
unpacked = {};    
ini = 1;
while ini < size(packed, 2)
    fin = find(isnan(packed(1, ini:end)), 1, 'first') + ini - 1;
    unpacked{end+1} = packed(:, ini:(fin-1))'; %#ok!
    ini = fin + 1;
    if ini >= size(packed, 2)
        break
    end
end


function packed = pack_polygons(unpacked)
% function packed = pack_polygons(unpacked)
% 
% Arguments:
%   unpacked = cell array or 2D matrices, each cell contains vertices for
%       one polygon with x-coord in row 1 and y-coord in row 2
% 
% Returns:
%   packed = polygon vertices as 2d array with points in columns and NaNs
%       separating each polygon
% % 
packed = [];
for ii = 1:length(unpacked)
   this_poly = unpacked{ii};
   packed = [packed, this_poly', nan(2,1)]; %#ok!
end