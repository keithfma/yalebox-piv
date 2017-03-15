function [xw, yw, xp, yp, done] = ...
    prep_world_coord_control_points(image_file, xw, yw, xp, yp, done)
% function [xw, yw, xp, yp, done] = ...
%     prep_world_coord_control_points(image_file, xw, yw, xp, yp, done)
%
% Create GUI to interactively define control points from the world coordinate
% image.
%
% NOTE: Providing the control points as input arguments initializes the data
%   table. This allows for restarting from a prior attempt or backup.
%
% NOTE: The table data (xw, yw, xp, yp, done) is backed up to backup.mat
%   after every GUI update. If something goes wrong, the GUI can be restarted
%   using this backup to recover.
%
% Arguments:
%   image_file: String, filename of the world coordinate grid image
%   xw, yw: Optional input, 1D vectors, control point world coordinate x- and y-
%       position, [m]
%   xp, yp: Optional input, 1D vectors, control point image coordinate x- and y-
%       position, [pixels]
%   done: Optional input, 1D vector, logical flag indicating whether point
%       definition is final (true) or not (false)
% 
% % Keith Ma

if nargin ~= 1 && nargin ~= 6
    error('Invalid number of input arguments');
end

%% constants

margin = 0.025; % norm
buffer = 0.01; % norm
large = 14; % pts
instruct = {'1. Enter world coordinates in xw, yw columns', ...
            '2. Select point in drop-down menu', ...
            '3. Push "Define" button to select pixel coordinates on image', ...
            '4. Mark as "done" if correct, otherwise redefine or delete', ...
            'NOTE: Points marked "done" cannot be edited or deleted', ...
            'NOTE: Progress is saved in a backup file to prevent disaster'};

im_left = 0.3 + margin;
im_bot = 0.1 + margin;
im_width = 1 - (margin + im_left);
im_height = 1 - (margin + im_bot);
im_pos = [im_left, im_bot, im_width, im_height];

tbl_left = margin;
tbl_bot = im_bot;
tbl_width = 1 - (2*margin + buffer + im_width);
tbl_height = 1 - (margin + tbl_bot);
tbl_pos = [tbl_left, tbl_bot, tbl_width, tbl_height];

list_left = margin;
list_bot = margin;
list_width = 0.1;
list_height = 0.05;
list_pos = [list_left, list_bot, list_width, list_height];

define_left = list_left + list_width + buffer;
define_bot = list_bot;
define_width = list_width;
define_height = list_height;
define_pos = [define_left, define_bot, define_width, define_height];

delete_left = define_left + define_width + buffer;
delete_bot = define_bot;
delete_width = define_width;
delete_height = define_height;
delete_pos = [delete_left, delete_bot, delete_width, delete_height];

done_left = delete_left + delete_width + buffer;
done_bot = delete_bot;
done_width = delete_width;
done_height = delete_height;
done_pos = [done_left, done_bot, done_width, done_height];

help_left = done_left + done_width + buffer;
help_bot = done_bot;
help_width = done_width;
help_height = done_height;
help_pos = [help_left, help_bot, help_width, help_height];

%% create GUI 

figure('Units', 'Normalized', 'Outerposition', [0 0 1 1], 'Tag', 'ctrl_gui');

axes('Units', 'Normalized', 'Position', im_pos, 'NextPlot', 'add', ...
    'Tag', 'ctrl_img');
imshow(imread(image_file));
title('World Coordinate Image');
hold on

spc = '               ';
uitable('Units', 'Normalized', ...
    'Position', tbl_pos, 'FontSize', large, ...
    'ColumnName', {['x [m]', spc], ['y [m]', spc], 'x [pix]', 'y [pix]', 'Done'}, ...
    'ColumnFormat', {'numeric', 'numeric', 'numeric', 'numeric', 'logical'}, ...
    'ColumnEditable', [true, true, false, false, true], ...
    'Tag', 'ctrl_table', 'CellEditCallback', @do_table);

uicontrol('Style', 'popupmenu', 'Units', 'normalized', ...
    'Position', list_pos, 'String', {'* SELECT *'}, 'FontSize', large, ...
    'Tag', 'ctrl_list');

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', define_pos, 'String', 'Define', 'FontSize', large, ...
    'Callback', @do_define);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', delete_pos, 'String', 'Delete', 'FontSize', large, ...
    'Callback', @do_delete);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', done_pos, 'String', 'Done', 'FontSize', large, ...
    'Callback', @do_done, 'Tag', 'ctrl_done');

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', help_pos, 'String', 'Help', 'FontSize', large, ...
    'Callback', @(~,~) msgbox(instruct));

% TODO

% set initial data (empty or from provided data)
if nargin == 1
    restore_state();
else
    restore_state(xw, yw, xp, yp, done);
end

update_all();

%% return results

% waitfor(findobj('Tag', 'ctrl_gui')); % waits here until the do_done() callback triggers uiresume
uiwait(findobj('Tag', 'ctrl_gui')); % waits here until the do_done() callback triggers uiresume
[data, ~] = get_table_data();
last = find(data(:, 5)==1, 1, 'last');
xw = data(1:last, 1);
yw = data(1:last, 2);
xp = data(1:last, 3);
yp = data(1:last, 4);
done = data(1:last, 5);
close(gcf);


function do_table(~, ~)

% identify edited cells
[curr, prev] = get_table_data;
[rr, cc] = find((curr ~= prev) & ~(isnan(curr) & isnan(prev)));

% revert illegal changes
if isscalar(cc) && isscalar(rr)
    % user edited one cell of xw, yw, done
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
    end
end

% update UI
set_table_data(curr);
update_all();

function [] = do_define(~, ~)
% Define image coordinates for selected world point, callback for define button
% % 

% get point index
idx = get_selected_point;
if isnan(idx); return; end

% interactive point selection
him = findobj('Tag', 'ctrl_img');
hpt = impoint(him);
pos = hpt.getPosition();
hpt.delete();
set(findobj('Tag', 'ctrl_gui'), 'WaitStatus', 'waiting'); % hack to preserve uiwait

% update table and plot
[data, ~] = get_table_data();
data(idx, 3:4) = pos;
set_table_data(data);
update_all();

function [] = do_delete(~, ~)
% Delete record for selected point, callback function for delete button
% %
idx = get_selected_point;
if isnan(idx); return; end

[curr, ~] = get_table_data();
curr(idx:end, :) = [curr(idx+1:end, :); nan(1, 4), false];
set_table_data(curr);
update_all();


function [] = do_done(~, ~)
% Callback for done button, deletes GUI which returns control to main function
% %

uiresume(findobj('Tag', 'ctrl_gui'));


function update_all()
% Update all GUI elements (point list, plot, done button) and save backup file
% %

update_point_list();
update_plot();
update_done_button();
backup_state();


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
if isnan(idx)
    warndlg('No point selected');
end

function update_plot()
% Update plotted markers on world coordinate image to match table
% %

% store view limits
him = findobj('Tag', 'ctrl_img');
xlim = him.XLim;
ylim = him.YLim;

% replot points and labels
delete(findobj('Type', 'text'));
delete(findobj('Type', 'line'));
[data, ~] = get_table_data();
for idx = 1:size(data, 1)
    xp = data(idx, 3);
    yp = data(idx, 4);
    if ~isnan(xp) && ~isnan(yp)
        text(xp, yp, num2str(idx), 'FontSize', 14, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
        plot(xp, yp, 'Color', 'r', 'Marker', 'x', 'MarkerSize', 12);
    end
end

% restore view limits
him.XLim = xlim;
him.YLim = ylim;


function update_done_button()
% Enable done button if ready to save results (i.e. no partial points)
% %

% label points
% % 0: no point, 1: partial, 2: complete not done, 3: complete done, 4: wtf
[data, ~] = get_table_data();
label = nan(size(data, 1), 1);
for idx = 1:size(data, 1)
    pts = data(idx, 1:4);
    done = data(idx, 5);
    if all(isnan(pts))
        label(idx) = 0;
    elseif any(isnan(pts)) && any(~isnan(pts))
        label(idx) = 1;
    elseif all(~isnan(pts)) && ~done
        label(idx) = 2;
    elseif all(~isnan(pts)) && done
        label(idx) = 3;
    else
        label(idx) = 4;
    end
end

% enable button if all points are complete or none
hd = findobj('Tag', 'ctrl_done');
if all(label == 0 | label == 3)
    hd.Enable = 'on';
else
    hd.Enable = 'off';
end


function backup_state()
% Save table data to backup file
% %

[data, ~] = get_table_data();
xw = data(:,1); %#ok!
yw = data(:,2); %#ok!
xp = data(:,3); %#ok!
yp = data(:,4); %#ok!
done = data(:,5); %#ok!
save('backup.mat', 'xw', 'yw', 'xp', 'yp', 'done');

function restore_state(xw, yw, xp, yp, done)
% Restore table data from backup file (if args are provided) or initialize
% new table (if args are not provided)
% %

if nargin == 0
    % init new table data
    num_pts = 500;
    data = [nan(num_pts, 4), false(num_pts, 1)];
else
    % restore table data
    num_pts = max(500, length(xw)+100); % make sure there is some extra space
    data = [nan(num_pts, 4), false(num_pts, 1)];
    data(1:length(xw), :) = [xw, yw, xp, yp, done];
end
set_table_data(data);
