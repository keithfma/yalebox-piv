function labels = prep_mask_labels(rgb, segments, labels)
% function labels = prep_mask_labels(rgb, segments, labels)
% 
% Interactive tool for defining super-pixel classes using ROI polygons. Select
% the region(s) that represent a single class (e.g. sand or other)
% 
% Arguments:
%   rgb: 3D, 24-bit color image
%   
%   segments: TODO
%
%   labels: TODO
% 
% Returns: 
%   labels: S
% % 

% NOTE: always interactive, we just save the training data, including the
% label vector.

% set defaults
narginchk(2, 3);
if nargin < 3; labels = []; end

% discard input labels if number of segments does not match
% note: this is a common case if the user updates segmentation settings
num_segments = length(unique(segments(:)));
if length(labels) ~= num_segments
    warning('Input labels do not match the number of image segments, ignoring them');
    labels = zeros(num_segments, 1, 'uint8');
end

% sanity checks: TODO

% TODO: we assume that segments starts with 1 and is sequential, check this

% constants
margin = 0.025; % norm
buffer = 0.01; % norm
font = 14; % pt
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

add_sand_left = margin;
add_sand_bot = margin;
add_sand_width = 0.1;
add_sand_height = 0.05;
add_sand_pos = [add_sand_left, add_sand_bot, add_sand_width, add_sand_height];

rem_sand_left = add_sand_left + add_sand_width + buffer;
rem_sand_bot = add_sand_bot;
rem_sand_width = add_sand_width;
rem_sand_height = add_sand_height;
rem_sand_pos = [rem_sand_left, rem_sand_bot, rem_sand_width, rem_sand_height];

add_other_left = rem_sand_left + rem_sand_width + buffer;
add_other_bot = rem_sand_bot;
add_other_width = rem_sand_width;
add_other_height = rem_sand_height;
add_other_pos = [add_other_left, add_other_bot, add_other_width, add_other_height];

rem_other_left = add_other_left + add_other_width + buffer;
rem_other_bot = add_other_bot;
rem_other_width = add_other_width;
rem_other_height = add_other_height;
rem_other_pos = [rem_other_left, rem_other_bot, rem_other_width, rem_other_height];

done_left = rem_other_left + rem_other_width + buffer;
done_bot = rem_other_bot;
done_width = rem_other_width;
done_height = rem_other_height;
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
    'Position', add_sand_pos, 'String', 'Add Sand', 'FontSize', font, ...
    'Callback', @do_add_sand);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Tag', 'AddBtn', ...
    'Position', rem_sand_pos, 'String', 'Remove Sand', 'FontSize', font, ...
    'Callback', @do_rem_sand);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Tag', 'NextBtn', ...
    'Position', add_other_pos, 'String', 'Add Other', 'FontSize', font, ...
    'Callback', @do_add_other);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Tag', 'NextBtn', ...
    'Position', rem_other_pos, 'String', 'Remove Other', 'FontSize', font, ...
    'Callback', @do_rem_other);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', done_pos, 'String', 'Done', 'FontSize', font, ...
    'Callback', @(~,~)uiresume);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', help_pos, 'String', 'Help', 'FontSize', font, ...
    'Callback', @(~,~)msgbox(instruct));

% display image
h_rgb = imshow(rgb);
set(h_rgb, 'Tag', 'RGB');
hold on

class_mask = labels(segments);
sand_mask = double(class_mask == 1);
other_mask = double(class_mask == 2);

[nr, nc, ~] = size(rgb);
blue = zeros(nr, nc, 3, 'uint8');
blue(:, :, 3) = 255;
h_sand = imshow(blue);
set(h_sand, 'AlphaData', 0.5*sand_mask, 'Tag', 'SandMask');

red = zeros(nr, nc, 3, 'uint8');
red(:, :, 1) = 255;
h_other = imshow(red); 
set(h_other, 'AlphaData', 0.5*other_mask, 'Tag', 'OtherMask');

hax.YDir = 'normal';

% store segments with figure for use in callbacks
setappdata(hf, 'segments', segments);
setappdata(hf, 'num_segments', length(unique(segments)));

% wait here until user presses "Done"
uiwait(); 

% gather results from transparency layers and cleanup
labels = zeros(num_segments, 1, 'uint8');

sand_mask = h_sand.AlphaData;
sand_idx = unique(segments(sand_mask > 0));
labels(sand_idx) = 1;

other_mask = h_other.AlphaData;
other_idx = unique(segments(other_mask > 0));
labels(other_idx) = 2;

close(hf)


function mask = select_by_polygon()
% Generate mask including all segments touched by user-drawn polygon
% %
num_segments = getappdata(gcf, 'num_segments');
segments = getappdata(gcf, 'segments');

poly = drawpolygon();
[nr, nc, ~] = size(get(findobj('Tag', 'RGB'), 'CData'));
footprint = poly2mask(poly.Position(:,1), poly.Position(:,2), nr, nc);
selected_segments = false(num_segments, 1);
selected_segments(unique(segments(footprint))) = true;
mask = selected_segments(segments);

delete(poly);


function [] = do_add_sand(~, ~)
% Callback for "Add Sand" button
% %
new_mask = select_by_polygon();
% add mask to sand
h_sand = findobj('Tag', 'SandMask');
sand_mask = h_sand.AlphaData;
sand_mask(new_mask) = 0.5;
h_sand.AlphaData = sand_mask;
% remove mask from other
h_other = findobj('Tag', 'OtherMask');
other_mask = h_other.AlphaData;
other_mask(new_mask) = 0;
h_other.AlphaData = other_mask;
% don't close window
uiwait();


function [] = do_rem_sand(~, ~)
% Callback for "Remove Sand" button
% %
new_mask = select_by_polygon();
% remove mask from sand
h_sand = findobj('Tag', 'SandMask');
sand_mask = h_sand.AlphaData;
sand_mask(new_mask) = 0;
h_sand.AlphaData = sand_mask;
% don't close window
uiwait();


function [] = do_add_other(~, ~)
% Callback for "Add Other" button
% %
new_mask = select_by_polygon();
% remove mask from sand
h_sand = findobj('Tag', 'SandMask');
sand_mask = h_sand.AlphaData;
sand_mask(new_mask) = 0;
h_sand.AlphaData = sand_mask;
% add mask to other
h_other = findobj('Tag', 'OtherMask');
other_mask = h_other.AlphaData;
other_mask(new_mask) = 0.5;
h_other.AlphaData = other_mask;
% don't close window
uiwait();

function [] = do_rem_other(~, ~)
% Callback for "Remove Other" button
% %
new_mask = select_by_polygon();
% remove mask from other
h_other = findobj('Tag', 'OtherMask');
other_mask = h_other.AlphaData;
other_mask(new_mask) = 0;
h_other.AlphaData = other_mask;
% don't close window
uiwait();
