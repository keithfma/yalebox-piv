function poly = prep_mask_define_labels(rgb)
% function poly = prep_mask_define_labels(rgb)
% 
% Interactive tool for defining pixel classes using ROI polygons. Select
% the region(s) that represent a single class (e.g. sand or other)
% 
% Arguments:
%   rgb: 3D, 24-bit color image
% 
% Returns: 
%  poly = polygon vertices as 2d array with points in columns and NaNs
%   separating each polygon
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

% create GUI, then wait until user presses "Done" button
hf = figure('Units', 'Normalized', 'Outerposition', [0 0 1 1]);

axes('Units', 'Normalized', 'Position', img_pos);

imshow(rgb);
set(gca, 'YDir', 'normal');

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', add_pos, 'String', 'Add', 'FontSize', font, ...
    'Callback', @do_add);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', done_pos, 'String', 'Done', 'FontSize', font, ...
    'Callback', @(~,~)uiresume);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', help_pos, 'String', 'Help', 'FontSize', font, ...
    'Callback', @(~,~)msgbox(instruct));

uiwait();

% retrieve, pack, and return polygons
rois = findobj(hf, 'Tag', 'ROI');
poly = [];
for ii = 1:length(rois)
   poly = [poly, rois(ii).Position', nan(2,1)]; %#ok!
end

close(hf)


function [] = do_add(~, ~)
% Callback for add mask layer button
% %
drawpolygon('Tag', 'ROI');
uiwait();
