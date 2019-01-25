function [] = prep_param()
% function [] = prep_param()
% 
% Create GUI for defining image prep parameters
% % 

update_path('util', 'prep')

%% constants

margin_left = 0.05;
margin_top = 0.025;

row_spc = 0.01;
row_height = 0.045;

col_spc = 0.0125;
col_width_1 = 0.25;
col_width_2 = 0.2;
col_width_3 = 0.3;
% col_width_4 = 0.15;
% col_width_5 = 0.25;
col_left_1 = margin_left;
col_left_2 = col_left_1 + col_width_1 + col_spc;
col_left_3 = col_left_2 + col_width_2 + col_spc;
% col_left_4 = col_left_3 + col_width_3 + col_spc;
% col_left_5 = col_left_4 + col_width_4 + col_spc;

%% create GUI
down_one = @(x) x - row_spc - row_height;

hf = figure('Units', 'Normalized', 'Outerposition', [0.25 0 0.5 1], 'Tag', 'prep_gui');
hf.Name = sprintf('Define Image Prep Parameters - Yalebox PIV - v%s', get_version());

param = prep_default_param;

row_bot = 1.0 - margin_top - row_height;
btn_param_file = uicontrol('Style', 'pushbutton', 'Units', 'Normalized');
btn_param_file.String = 'Parameter File';
btn_param_file.Position = [col_left_1, row_bot, col_width_1, row_height];
btn_param_file.TooltipString = 'select new or existing parameter file';
btn_param_file.Callback = @btn_param_file_callback;

txt_param_file = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
txt_param_file.Position = [col_left_2, row_bot, col_width_2 + col_width_3, row_height];
txt_param_file.Tag = 'txt_param_file';

row_bot = down_one(row_bot);  % new row -----------------------------------
 
btn_image_dir = uicontrol('Style', 'pushbutton', 'Units', 'Normalized');
btn_image_dir.String = 'Image Directory';
btn_image_dir.Position = [col_left_1, row_bot, col_width_1, row_height];
btn_image_dir.TooltipString = param.images.path.help;
btn_image_dir.Callback = @btn_image_dir_callback;

txt_image_dir = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
txt_image_dir.Position = [col_left_2, row_bot, col_width_2 + col_width_3, row_height];
txt_image_dir.Tag = 'txt_image_dir';

row_bot = down_one(row_bot);  % new row -----------------------------------

btn_woco_file = uicontrol('Style', 'pushbutton', 'Units', 'Normalized');
btn_woco_file.String = 'World Coordinate File';
btn_woco_file.Position = [col_left_1, row_bot, col_width_1, row_height];
btn_woco_file.TooltipString = param.images.woco_file.help;

row_bot = down_one(row_bot);  % new row -----------------------------------

btn_test_file = uicontrol('Style', 'pushbutton', 'Units', 'Normalized');
btn_test_file.String = 'Test Image File';
btn_test_file.Position = [col_left_1, row_bot, col_width_1, row_height];
btn_test_file.TooltipString = param.test.test_file.help;

row_bot = down_one(row_bot);  % new row -----------------------------------

btn_def_woco = uicontrol('Style', 'pushbutton', 'Units', 'Normalized');
btn_def_woco.String = 'World Coordinate Control Points';
btn_def_woco.Position = [col_left_1, row_bot, col_width_1, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

btn_run_crop = uicontrol('Style', 'pushbutton', 'Units', 'Normalized');
btn_run_crop.String = 'Rectify and Crop';
btn_run_crop.Position = [col_left_1, row_bot, col_width_1, row_height];

txt_crop_xlim = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_crop_xlim.String = 'X Min, Max (m)';
txt_crop_xlim.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_crop_xlim = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_crop_xlim.String = '';
edit_crop_xlim.Position = [col_left_3, row_bot, col_width_3, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

txt_crop_ylim = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_crop_ylim.String = 'Y Min, Max (m)';
txt_crop_ylim.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_crop_ylim = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_crop_ylim.String = '';
edit_crop_ylim.Position = [col_left_3, row_bot, col_width_3, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

btn_def_mask_manual = uicontrol('Style', 'pushbutton', 'Units', 'Normalized');
btn_def_mask_manual.String = 'Manual Mask';
btn_def_mask_manual.Position = [col_left_1, row_bot, col_width_1, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

btn_run_auto_mask = uicontrol('Style', 'pushbutton', 'Units', 'Normalized');
btn_run_auto_mask.String = 'Automatic Mask';
btn_run_auto_mask.Position = [col_left_1, row_bot, col_width_1, row_height];

txt_mask_hue_lim = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_mask_hue_lim.String = 'Hue Min, Max';
txt_mask_hue_lim.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_mask_hue_lim = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_mask_hue_lim.String = '';
edit_mask_hue_lim.Position = [col_left_3, row_bot, col_width_3, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

txt_mask_val_lim = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_mask_val_lim.String = 'Value Min, Max';
txt_mask_val_lim.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_mask_val_lim = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_mask_val_lim.String = '';
edit_mask_val_lim.Position = [col_left_3, row_bot, col_width_3, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

txt_mask_ent_lim = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_mask_ent_lim.String = 'Entropy Min, Max';
txt_mask_ent_lim.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_mask_ent_lim = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_mask_ent_lim.String = '';
edit_mask_ent_lim.Position = [col_left_3, row_bot, col_width_3, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

txt_mask_ent_len = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_mask_ent_len.String = 'Entropy Window Size (pixels)';
txt_mask_ent_len.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_mask_ent_len = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_mask_ent_len.String = '';
edit_mask_ent_len.Position = [col_left_3, row_bot, col_width_3, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

txt_mask_morph_open = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_mask_morph_open.String = 'Morphological Opening Radius (pixels)';
txt_mask_morph_open.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_mask_morph_open = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_mask_morph_open.String = '';
edit_mask_morph_open.Position = [col_left_3, row_bot, col_width_3, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

txt_mask_morph_erode = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_mask_morph_erode.String = 'Morphological Erosion Radius (pixels)';
txt_mask_morph_erode.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_mask_morph_erode = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_mask_morph_erode.String = '';
edit_mask_morph_erode.Position = [col_left_3, row_bot, col_width_3, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

btn_run_hist = uicontrol('Style', 'pushbutton', 'Units', 'Normalized');
btn_run_hist.String = 'Histogram Equalization';
btn_run_hist.Position = [col_left_1, row_bot, col_width_1, row_height];

txt_eql_len = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_eql_len.String = 'Equalization Window Size (pixels)';
txt_eql_len.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_eql_len = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_eql_len.String = '';
edit_eql_len.Position = [col_left_3, row_bot, col_width_3, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

btn_run_series = uicontrol('Style', 'pushbutton', 'Units', 'Normalized');
btn_run_series.String = 'Test Series';
btn_run_series.Position = [col_left_1, row_bot, col_width_1, row_height];

txt_test_num_img = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_test_num_img.String = 'Num Images';
txt_test_num_img.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_test_num_img = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_test_num_img.String = '';
edit_test_num_img.Position = [col_left_3, row_bot, col_width_3, row_height];

row_bot = down_one(row_bot);  % new row -----------------------------------

txt_test_out_file = uicontrol('Style', 'text', 'Units', 'Normalized', 'HorizontalAlignment', 'right');
txt_test_out_file.String = 'Output File';
txt_test_out_file.Position = [col_left_2, row_bot, col_width_2, row_height];

edit_test_out_file = uicontrol('Style', 'edit', 'Units', 'Normalized', 'HorizontalAlignment', 'left');
edit_test_out_file.String = '';
edit_test_out_file.Position = [col_left_3, row_bot, col_width_3, row_height];

% store default parameters and update the various fields
guidata(hf, param);
update_view(hf);

end


function update_file(hobj)
% (Re)save parameters to specified file
% %

param = guidata(hobj);
txt_param_file = findobj('Tag', 'txt_param_file');
param_file = txt_param_file.String;
if ~isempty(param_file)
    save(param_file, 'param');
end
% TODO: else alert missing
end


function update_view(hobj)
% Update displayed fields using parameters stored in guidata
% %
param = guidata(hobj);

txt_image_dir = findobj('Tag', 'txt_image_dir');
txt_image_dir.String = param.images.path.value;

end


function btn_param_file_callback(hobject, eventdata)
% Set parameter file name
% %

filename = uiputfile('*.mat');
txt_param_file = findobj('Tag', 'txt_param_file');
txt_param_file.String = filename;
    
end


function btn_image_dir_callback(hobject, eventdata)
% Set image directory
% %

param = guidata(hobject);
path = uigetdir('Select directory containing image files');
path = strip(path, 'right', filesep);
param.images.path.value = path;
guidata(hobject, param);
update_view(hobject);
update_file(hobject);

end

