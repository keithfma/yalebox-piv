function [] = new_ctrl_pts()
% new, simpler program for collecting world-coordinate control points

% The main function creates the GUI and registers all the callbacks


% debug: cleanup and define fake inputs
close all
image_file = '..\data\fault_ss_03_siden_clean\woco\fault_ss_03_siden_woco_b.JPG';

% constants
margin = 0.05;
large = 14;
buffer = 0.1;

im_left = 0.3 + margin;
im_bot = 0.1 + margin;
im_width = 1 - (margin + im_left);
im_height = 1 - (margin + im_bot);
im_pos = [im_left, im_bot, im_width, im_height];

tbl_left = margin;
tbl_bot = im_bot;
tbl_width = 1 - (3*margin + im_width);
tbl_height = 1 - (margin + tbl_bot);
tbl_pos = [tbl_left, tbl_bot, tbl_width, tbl_height];

text_left = margin;
text_bot = margin;
text_width = 0.1;
text_height = 0.05;
text_pos = [text_left, text_bot, text_width, text_height];

curr_left = text_left + text_width*(1 + buffer);
curr_bot = text_bot;
curr_width = text_width;
curr_height = text_height;
curr_pos = [curr_left, curr_bot, curr_width, curr_height];

define_left = curr_left + (1 + buffer)*curr_width;
define_bot = curr_bot;
define_width = curr_width;
define_height = curr_height;
define_pos = [define_left, define_bot, define_width, define_height];

delete_left = define_left + 1.1*define_width;
delete_bot = define_bot;
delete_width = define_width;
delete_height = define_height;
delete_pos = [delete_left, delete_bot, delete_width, delete_height];

done_left = delete_left + 1.1*delete_width;
done_bot = delete_bot;
done_width = delete_width;
done_height = delete_height;
done_pos = [done_left, done_bot, done_width, done_height];

% create GUI
figure('Units', 'Normalized', 'Outerposition', [0 0 1 1]);

axes('Units', 'Normalized', 'Position', im_pos, ...
    'Tag', 'ctrl_img');
imshow(imread(image_file));
title('World Coordinate Image');

uitable('Data', cell(100, 5), 'Units', 'Normalized', ...
    'Position', tbl_pos, ...
    'ColumnName', {'x_world', 'y_world', 'x_pixel', 'y_pixel', 'done'}, ...
    'ColumnFormat', {'numeric', 'numeric', 'numeric', 'numeric', 'logical'}, ...
    'ColumnEditable', [true, true, false, false, true], ...
    'Tag', 'ctrl_table', 'CellEditCallback', @update_point_list);

uicontrol('Style', 'text', 'Units', 'normalized', ...
    'Position', text_pos, 'String', 'Control Point:', 'FontSize', large);

uicontrol('Style', 'popupmenu', 'Units', 'normalized', ...
    'Position', curr_pos, 'String', {'Current', 'Next'}, 'FontSize', large);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', define_pos, 'String', 'Define', 'FontSize', large, ...
    'Callback', @do_define);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', delete_pos, 'String', 'Delete', 'FontSize', large, ...
    'Callback', @do_delete);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', done_pos, 'String', 'Done', 'FontSize', large, ...
    'Callback', @do_done);

function [] = do_define(~, ~)
% Callback function for next button
warndlg('DO DEFINE');

function [] = do_delete(~, ~)
% Callback function for prev button
warndlg('DO DELETE');

function [] = do_done(~, ~)
% Callback for done button
warndlg('DO DONE');

function update_point_list(ht, ~)
% Populate point list with definable points, callback for table edit
warndlg('DO UPDATE POINT LIST');   