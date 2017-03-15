function mask = prep_mask_manual_new(img)
% function mask = prep_mask_manual(img)
%
% Interactive GUI to define a mask to black out region(s) of an image.
%
% Arguments:
%   img = Matrix, image data that can be displayed by imshow
%   mask = 2D matrix, logical, true where there is sand and false elsewhere.
%
% % Keith Ma

%% constants

margin = 0.025; % norm
buffer = 0.01; % norm
font = 14; % pt
instruct = {'Some', 'Instructions', 'Here'};

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

data = struct('img', img, 'mask', true(size(img, 1), size(img, 2)), 'layer', {{}});
axes('Units', 'Normalized', 'Position', img_pos, 'NextPlot', 'add', ...
    'Tag', 'mask_img', 'UserData', data);
image(img);
axis equal off

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

uiwait(findobj('Tag', 'mask_gui'));

%% return results

hi = findobj('Tag', 'mask_img'); 
data = hi.UserData;
mask = data.mask;
close(findobj('Tag', 'mask_gui'));


function do_add(~, ~)
% Callback for add mask layer button
% %

% get current data
hi = findobj('Tag', 'mask_img');
data = hi.UserData;

% get new mask polygon interactively
p = impoly('Closed', true);
set(findobj('Tag', 'mask_gui'), 'WaitStatus', 'waiting'); % hack for uiwait
this_layer = ~createMask(p); % invert
data.layer{end+1} = this_layer;

% compose mask 
mask = true(size(data.mask));
for ii = 1:length(data.layer)
    mask = mask & data.layer{ii};
end
data.mask = mask;

% updated stored data
hi.UserData = data;

% update plot
masked_img = data.img;
masked_img(repmat(~mask, 1, 1, 3)) = 0;
axes(hi);
delete(findobj('Type', 'Image', 'Parent', hi));
image(masked_img);
axis equal off


% TODO: add delete button

function do_done(~, ~)
% Callback for done button, disables uiwait
% %
uiresume(findobj('Tag', 'mask_gui'));