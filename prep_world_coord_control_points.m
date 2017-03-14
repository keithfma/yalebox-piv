function [] = prep_world_coord_control_points(image_file, output_file, backup_file)
% function [] = prep_world_coord_control_pts(image_file, output_file, backup_file)
%
% Create GUI to interactively define control points from the world coordinate
% image. Results are saved to a .mat file containing four variables: ctrl_xw,
% ctrl_yw, ctrl_xp, ctrl_yp (defined below).
%
% NOTE: Results are saved to .mat file after every user update (see the
% save_state function for the backup file name). If something goes wrong, the
% GUI can be restarted using this backup to recover.
%
% NOTE: This GUI replaces a prior version that used a different syntax. If you
% see crazy errors, then the calling program may be trying to use the old
% version.
%
% Arguments:
%   image_file: String, filename of the world coordinate grid image
%   output_file: String, filename of .mat file to save results, note that this
%       value can be changed in the GUI
%   backup_file: [Optional] String, filename of previous results or backup to be
%       loaded initially, this allows for restarting if something goes wrong.
%   ctrl_xw, ctrl_yw: 1D vectors, control point world coordinate x- and
%       y-position in meters
%   ctrl_xp, ctrl_yp = 1D vectors, control point image coordinate x- and y-position
%       in pixels
% 
% % Keith Ma

% constants --------------------------------------------------------------------

margin = 0.025;
large = 14;
buffer = 0.1;
max_pts = 100;

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

delete_left = define_left + (1 + buffer)*define_width;
delete_bot = define_bot;
delete_width = define_width;
delete_height = define_height;
delete_pos = [delete_left, delete_bot, delete_width, delete_height];

done_left = delete_left + (1 + buffer)*delete_width;
done_bot = delete_bot;
done_width = delete_width;
done_height = delete_height;
done_pos = [done_left, done_bot, done_width, done_height];

name_left = done_left + (1 + buffer)*done_width;
name_bot = done_bot;
name_width = done_width*2;
name_height = done_height;
name_pos = [name_left, name_bot, name_width, name_height];

% create GUI -------------------------------------------------------------------

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

uicontrol('Style', 'edit', 'Units', 'normalized', ...
    'Position', name_pos, 'String', output_file, 'FontSize', large, ...
    'Tag', 'ctrl_file');

% set initial data (empty or from backup file)
if nargin < 3 || isempty(backup_file)
    init_data = [nan(max_pts, 4), false(max_pts, 1)];
else
    load(backup_file, 'ctrl_xw', 'ctrl_yw', 'ctrl_xp', 'ctrl_yp');
    max_pts = length(ctrl_xw);
    init_data = [ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, false(max_pts, 1)];
end
set_table_data(init_data);

update_all();

% define callbacks ------------------------------------------------------------

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
% Callback for done button, check points and write results
% %

hf = findobj('Tag', 'ctrl_file');
filename = hf.String;
if ~isempty(filename)
    save_state(filename)
    close(findobj('Tag', 'ctrl_gui'));
else
    warndlg('Cannot save results without an output file name');
end

% define utilities -------------------------------------------------------------

function update_all()
% Update all GUI elements (point list, plot, done button) and save backup file
% %

update_point_list();
update_plot();
update_done_button();
save_state();


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
hf = findobj('Tag', 'ctrl_file');
if all(label == 0 | label == 3)
    hd.Enable = 'on';
    hf.Enable = 'on';
else
    hd.Enable = 'off';
    hf.Enable = 'off';
end

function save_state(filename)
% Save the current data table, saves backup if filename is empty, deletes backup
% on final save
% %

backup_file = 'ctrl_pts_backup.mat';

if nargin == 0 || isempty(filename)
    filename = backup_file;
else
    delete(backup_file);
end

[data, ~] = get_table_data();
ctrl_xw = data(:, 1); %#ok!
ctrl_yw = data(:, 2); %#ok!
ctrl_xp = data(:, 3); %#ok!
ctrl_yp = data(:, 4); %#ok!
save(filename, 'ctrl_xw', 'ctrl_yw', 'ctrl_xp', 'ctrl_yp');