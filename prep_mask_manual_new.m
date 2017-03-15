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

axes('Units', 'Normalized', 'Position', img_pos, 'NextPlot', 'add', ...
    'Tag', 'mask_img');
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

%% return results


