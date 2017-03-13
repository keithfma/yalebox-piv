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
max_pts = 5;

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

list_left = margin;
list_bot = margin;
list_width = 0.1;
list_height = 0.05;
list_pos = [list_left, list_bot, list_width, list_height];

define_left = list_left + (1 + buffer)*list_width;
define_bot = list_bot;
define_width = list_width;
define_height = list_height;
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

uitable('Units', 'Normalized', ...
    'Position', tbl_pos, 'FontSize', large, ...
    'ColumnName', {'x [m]', 'y [m]', 'x [pix]', 'y [pix]', 'Done'}, ...
    'ColumnFormat', {'numeric', 'numeric', 'numeric', 'numeric', 'logical'}, ...
    'ColumnEditable', [true, true, false, false, true], ...
    'Tag', 'ctrl_table', 'CellEditCallback', @do_update_table);

uicontrol('Style', 'popupmenu', 'Units', 'normalized', ...
    'Position', list_pos, 'String', {'-SELECT-'}, 'FontSize', large, ...
    'Tag', 'ctrl_list');

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', define_pos, 'String', 'Define', 'FontSize', large, ...
    'Callback', @do_define);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', delete_pos, 'String', 'Delete', 'FontSize', large, ...
    'Callback', @do_delete);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', done_pos, 'String', 'Done', 'FontSize', large, ...
    'Callback', @do_done);

set_table_data([[1, 2, 3, 4, false]; nan(max_pts, 4), false(max_pts, 1)]);
update_point_list();


function set_table_data(data)
% Set current and update previous data in control point table
% 
% Converts data matrix to cell array, preserving logical class, and updates both
% the current data (table.Data) and previous data (table.UserData)
% %
data_cell = [num2cell(data(:, 1:4)), num2cell(logical(data(:, 5)))];
ht = findobj('Tag', 'ctrl_table');
ht.Data = data_cell;
ht.UserData = data_cell;


function [curr, prev] = get_table_data()
% Fetch current and previous data from control point table
% 
% Converts stored cell arrays to matrixes, converting all elements to double
% %
ht = findobj('Tag', 'ctrl_table');
curr = [cell2mat(ht.Data(:,1:4)), cell2mat(ht.Data(:, 5))];
prev = [cell2mat(ht.UserData(:,1:4)), cell2mat(ht.UserData(:, 5))];


function update_point_list()
% Populate point list with definable points (have xw, yw, and not marked "done")
% %
hl = findobj('Tag', 'ctrl_list');
active_points = hl.String(1); % keep prompt
[data, ~] = get_table_data();
for ii = 1:size(data, 1)
    if ~isnan(data(ii,1)) && ~isnan(data(ii,2)) && ~data(ii,5) 
        active_points{end+1} = ii; %#ok!
    end
end
hl.String = active_points;
hl.Value = length(active_points);


function idx = get_selected_point()
% Return index of point selected in point list, or NaN if none
% %
hl = findobj('Tag', 'ctrl_list');
idx = str2double(hl.String{hl.Value});


function do_update_table(~, ~)
% Callback for table update

% fetch current and previous table state as matrix
[curr, prev] = get_table_data;

% identify edited cell
[rr, cc] = find((curr ~= prev) & ~(isnan(curr) & isnan(prev)));

% validate user changes
if cc == 5
    % update "done" state
    if any(isnan(curr(rr, 1:4)))
        warndlg('Cannot mark incomplete control point as "done"');
        curr(rr, cc) = prev(rr, cc); % undo
    end
elseif cc == 1 || cc == 2
    % update world coordinates
    if curr(rr, 5)
        warndlg('Cannot update world coordinate for point marked as "done"');
        curr(rr, cc) = prev(rr, cc);  % undo
    end
else
    % this should not be possible
    warndlg('This should not be possible, undoing last table update');
    curr(rr, cc) = prev(rr, cc);  % undo
end

% update UI
set_table_data(curr);
update_point_list();


function [] = do_define(~, ~)
% Define image coordinates for selected world point, callback for define button
% % 
idx = get_selected_point;
if isnan(idx)
    warndlg('Cannot define: no point selected');
    return
end
% TODO: interactive point selection
warndlg(sprintf('DEFINE: %d', idx));


function [] = do_delete(~, ~)
% Delete record for selected point, callback function for delete button
idx = get_selected_point;
if isnan(idx)
    warndlg('Cannot delete: No point selected');
    return
end
[curr, ~] = get_table_data();
curr(idx:end, :) = [curr(idx+1:end, :); nan(1, 4), false];
set_table_data(curr);
update_point_list();


function [] = do_done(~, ~)
% Callback for done button
% TODO: confirm before exiting
warndlg('DO DONE');