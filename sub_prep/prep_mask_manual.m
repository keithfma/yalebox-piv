function polygons = prep_mask_manual(img, polygons)
% function polygons = prep_mask_manual(img, polygons) % TODO: accept existing
% polygons
%
% Interactive GUI to define a mask to black out region(s) of an image.
%
% Arguments:
%   img = Matrix, image data that can be displayed by imshow polygons =
%   polygons: optional 2D array, vertices of mask polygons, x-coords in row 1 and
%       y-coords in row 2, polygons separated by NaN
%
% % Keith Ma

% TODO: also store and return the mask, needed for series processing

% set defaults
if nargin < 2; polygons = []; end

% sanity check
narginchk(1,2);
validateattributes(img, {'numeric'}, {'3d'});
validateattributes(polygons, {'numeric'}, {'2d'});

%% constants

margin = 0.025; % norm
buffer = 0.01; % norm
font = 14; % pt
instruct = {'1. Click "Add" button to define a mask polygon interactively', ...
            '2. Click "Delete" button to remove the most recent mask layer', ...
            '3. Click "Done" to return'};

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

del_left = add_left + add_width + buffer;
del_bot = add_bot;
del_width = add_width;
del_height = add_height;
del_pos = [del_left, del_bot, del_width, del_height];

done_left = del_left + del_width + buffer;
done_bot = del_bot;
done_width = del_width;
done_height = del_height;
done_pos = [done_left, done_bot, done_width, done_height];

help_left = done_left + done_width + buffer;
help_bot = done_bot;
help_width = done_width;
help_height = done_height;
help_pos = [help_left, help_bot, help_width, help_height];


%% create GUI

figure('Units', 'Normalized', 'Outerposition', [0 0 1 1], 'Tag', 'mask_gui');

data = struct(...
    'img', img, ...
    'polygons', {unpack_polygons(polygons)});

axes('Units', 'Normalized', 'Position', img_pos, 'NextPlot', 'add', ...
    'Tag', 'mask_img', 'UserData', data);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', add_pos, 'String', 'Add New Mask', 'FontSize', font, ...
    'Callback', @do_add);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', del_pos, 'String', 'Delete Last Mask', 'FontSize', font, ...
    'Callback', @do_del);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', done_pos, 'String', 'Done', 'FontSize', font, ...
    'Callback', @do_done);

uicontrol('Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', help_pos, 'String', 'Help', 'FontSize', font, ...
    'Callback', @(~,~)msgbox(instruct));

update_mask();
uiwait(findobj('Tag', 'mask_gui'));

%% return results

hi = findobj('Tag', 'mask_img'); 
data = hi.UserData;
polygons = pack_polygons(data.polygons);
close(findobj('Tag', 'mask_gui'));


function packed = pack_polygons(unpacked)
% Pack polygon vertices as 2d array with points in columns and NaNs
% separating each polygon
% % 
packed = [];
for ii = 1:length(unpacked)
   this_poly = unpacked{ii};
   packed = [packed, this_poly', nan(2,1)]; %#ok!
end


function unpacked = unpack_polygons(packed)
% Pack polygon vertices as 2d array with points in columns and NaNs
% separating each polygon
% % 
unpacked = {};    
ini = 1;
while ini < size(packed, 2)
    fin = find(isnan(packed(1, ini:end)), 1, 'first') + ini - 1;
    disp([packed(1,ini), packed(1,fin)]);
    unpacked{end+1} = packed(:, ini:(fin-1))'; %#ok!
    ini = fin + 1;
    if ini >= size(packed, 2)
        break
    end
end


function do_add(~, ~)
% Callback for add mask layer button
% %

% get new mask polygon interactively
this_poly = impoly(gca, 'Closed', true);
set(findobj('Tag', 'mask_gui'), 'WaitStatus', 'waiting'); % hack for uiwait
if isempty(this_poly)
    return % aborted
end

% store data and update GUI
hi = findobj('Tag', 'mask_img');
data = hi.UserData;
data.polygons{end+1} = getPosition(this_poly);
hi.UserData = data;
delete(this_poly);
update_mask();


function update_mask()
% Refresh mask and plot
% %

% fetch current data
hi = findobj('Tag', 'mask_img');
data = hi.UserData;

% ensure the image is displayed so mask can be created
if isempty(hi.Children)
    image(data.img);
end

% compose mask 
mask = true(size(data.img,1), size(data.img, 2));
for ii = 1:length(data.polygons)
    this_poly = impoly(gca, data.polygons{ii}, 'Closed', true);
    mask = mask & ~createMask(this_poly);
    delete(this_poly);
end

% replot masked image
delete(hi.Children);
masked_img = data.img;
masked_img(repmat(~mask, 1, 1, 3)) = 0;
axes(hi);
delete(findobj('Type', 'Image', 'Parent', hi));
image(masked_img);
axis equal off

hi.UserData = data;


function do_del(~, ~)
% Callback function for delete button
% %
hi = findobj('Tag', 'mask_img');
data = hi.UserData;
if ~isempty(data.polygons)
    data.polygons = data.polygons(1:end-1);    
end
hi.UserData = data;
update_mask();


function do_done(~, ~)
% Callback for done button, disables uiwait
% %
uiresume(findobj('Tag', 'mask_gui'));