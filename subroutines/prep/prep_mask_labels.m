function [labels, poly_out] = prep_mask_labels(rgb, poly_in, interactive)
% function [labels, polygons] = prep_mask_labels(rgb, poly_in, interactive)
% 
% Interactive tool for defining pixel classes using ROI polygons. Select
% the region(s) that represent a single class (e.g. sand or other)
% 
% Arguments:
%   rgb: 3D, 24-bit color image
%
%   poly_in: initial ROI polygons for labeled image area, cell array of
%       2D matrices, each cell contains vertices for one polygon with
%       x-coord in column 1 and y-coord in column 2
%
%   interactive: logical scalar, set true to define/update label polygons
%       interactively, default is false
% 
% Returns: 
%   labels: integer array, linear indices into one layer of RGB for all
%       pixels inside the label polygons
% 
%   poly_out = final ROI polygons for labeled image area,  cell array of
%       2D matrices, each cell contains vertices for one polygon with
%       x-coord in column 1 and y-coord in column 2
% % 

% set defaults
if nargin < 3; interactive = false; end

% TODO: sanity checks

if interactive
    % interactive polygon definition
    poly_out = define_polygons(rgb, poly_in);
else
    % use input polygons without change
    poly_out = poly_in;
end

% TODO: do multiclass labeling
mask = false(size(rgb, 1), size(rgb, 2));
for ii = 1:length(poly_out)
    x_pts = poly_out{ii}(:, 1);
    y_pts = poly_out{ii}(:, 2);
    mask = mask | poly2mask(x_pts, y_pts, size(mask, 1), size(mask, 2));
end 
labels = find(mask);


function outpoly = define_polygons(rgb, inpoly)
% Interactive GUI for defining/updating label ROI polygons
% 
% Arguments:
%   rgb: 3D, 24-bit color image
% 
%   inpoly = cell array of 2D matrices, each cell contains vertices for one
%       polygon with x-coord in column 1 and y-coord in column 2
% 
% Returns:
%   outpoly = cell array of 2D matrices, each cell contains vertices for
%       one polygon with x-coord in column 1 and y-coord in column 2
% %

% constants
margin = 0.025; % norm
buffer = 0.01; % norm
font = 14; % pt
% TODO: update instrictions
instruct = {'Create polygons for each training class - the first class is sand', ...
            '+ To create a polygon: click "Add" button and then click to define vertices', ...
            '+ To change classes: click "Prev" or "Next" buttons', ...
            '+ To edit a polygon: click and vertices', ...
            '+ To delete a polygon: right click and select "Delete Polygon"', ...
            '+ Finish: click "Done"'};

img_left = margin;
img_bot = 0.1;
img_width = 1 - 2*margin;
img_height = 1 - (img_bot + 2*margin);
img_pos = [img_left, img_bot, img_width, img_height];

prev_left = margin;
prev_bot = margin;
prev_width = 0.1;
prev_height = 0.05;
prev_pos = [prev_left, prev_bot, prev_width, prev_height];

add_left = prev_left + prev_width + buffer;
add_bot = prev_bot;
add_width = prev_width;
add_height = prev_height;
add_pos = [add_left, add_bot, add_width, add_height];

next_left = add_left + add_width + buffer;
next_bot = add_bot;
next_width = add_width;
next_height = add_height;
next_pos = [next_left, next_bot, next_width, next_height];

done_left = next_left + next_width + buffer;
done_bot = next_bot;
done_width = next_width;
done_height = next_height;
done_pos = [done_left, done_bot, done_width, done_height];

help_left = done_left + done_width + buffer;
help_bot = done_bot;
help_width = done_width;
help_height = done_height;
help_pos = [help_left, help_bot, help_width, help_height];

% create UI
hf = figure('Units', 'Normalized', 'Outerposition', [0 0 1 1]);

hax = axes('Units', 'Normalized', 'Position', img_pos);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Tag', 'PrevBtn', ...
    'Position', prev_pos, 'String', 'Prev Class', 'FontSize', font, ...
    'Callback', @do_prev);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Tag', 'AddBtn', ...
    'Position', add_pos, 'String', 'Add Poly', 'FontSize', font, ...
    'Callback', @do_add);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Tag', 'NextBtn', ...
    'Position', next_pos, 'String', 'Next Class', 'FontSize', font, ...
    'Callback', @do_next);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', done_pos, 'String', 'Done', 'FontSize', font, ...
    'Callback', @(~,~)uiresume);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', help_pos, 'String', 'Help', 'FontSize', font, ...
    'Callback', @(~,~)msgbox(instruct));

colors = lines;
setappdata(hf, 'class_colors', colors);
setappdata(hf, 'prev_class', 1);
setappdata(hf, 'current_class', 1);
setappdata(hf, 'next_class', 2);
update_colors();

% display image and input polygons 
imshow(rgb);
for ii = 1:length(inpoly)
    images.roi.Polygon(gca, ...
        'Tag', 'Class', ...
        'Position', inpoly{ii}, ...
        'Color', colors(ii, :), ...
        'UserData', ii);
end
hax.YDir = 'normal';

% wait here until user presses "Done"
uiwait(); 

% % retreive and reformat output polygons
objs = findobj(hf, 'Tag', 'ROI');

class_ids = [];
for jj = 1:length(objs)
    class_ids(end+1) = objs(jj).UserData;  %#ok<AGROW>
end
    
outpoly = cell(size(unique(class_ids)));
for jj = 1:length(outpoly)
    outpoly{jj} = {};
end

for kk = 1:length(objs)
    class_id = objs(kk).UserData;
    outpoly{class_id}{end+1} = objs(kk).Position;
end

close(hf)


function [] = update_colors()
% Utility for updating colors when current class changes
% %
colors = getappdata(gcf, 'class_colors');

prev_obj = findobj(gcf, 'Tag', 'PrevBtn');
prev_class = getappdata(gcf, 'prev_class');
prev_obj.ForegroundColor = colors(prev_class, :);

add_obj = findobj(gcf, 'Tag', 'AddBtn');
current_class = getappdata(gcf, 'current_class');
add_obj.ForegroundColor = colors(current_class, :);

next_obj = findobj(gcf, 'Tag', 'NextBtn');
next_class = getappdata(gcf, 'next_class');
next_obj.ForegroundColor = colors(next_class, :);


function [] = update_class(delta)
% Update class indices
% %
current_class = getappdata(gcf, 'current_class');
current_class = max(current_class + delta, 1);
prev_class = max(current_class - 1, 1);
next_class = current_class + 1;
setappdata(gcf, 'prev_class', prev_class);
setappdata(gcf, 'current_class', current_class);
setappdata(gcf, 'next_class', next_class);
update_colors();


function [] = do_prev(~, ~)
% Callback for "Prev Class" button
% %
update_class(-1);

function [] = do_add(~, ~)
% Callback for "Add Poly" button
% %
colors = getappdata(gcf, 'class_colors');
current_class = getappdata(gcf, 'current_class');
drawpolygon('Tag', 'ROI', 'Color', colors(current_class, :), 'UserData', current_class);
uiwait();


function [] = do_next(~, ~)
% Callback for "Next Class" button
% %
update_class(1)
