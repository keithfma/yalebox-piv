% Template script for exploring and selecting image pre-processing parameters.
% The intended wusage is to make a copy of the script for a given experiment,
% run it cell by cell, modifying the default parameters to suit the experiment
% particulars.
%
% NOTE: Specific values are from experiment fault_ss_03_sidef, they may be
% useful or not depending on the experiment

%% Init

% define
image_path = '../data/fault_ss_03_sidef_clean'; % linux slash always works
image_glob = [image_path, '/fault_ss_03_sidef_*.jpg'];
woco_image_file = [image_path, '/woco/fault_ss_03_sidef_woco_b.JPG'];
step_image_file = [image_path, '/fault_ss_03_sidef_613.jpg'];
num_test_images = 5;
test_result_file = 'junk.nc'; % netCDF
ctrl_pts_file = 'fault_ss_03_sidef_prep_ctrl_pts.mat'; % MAT
param_file =  'fault_ss_03_sidef_prep_param.mat'; % MAT
result_file = '../data/fault_ss_03_sidef_image.nc'; % netCDF

% get image file name list
tmp = dir(image_glob);
image_names = {tmp(:).name};

% pick well-distributed random test images
edge = round(linspace(1, length(image_names), num_test_images+1));
start = edge(1:end-1);
stop = edge(2:end);
idx = round(rand(size(start)).*(stop - start) + start);
test_image_names = image_names(idx);

%% Define coordinate system control points

mode = 'load'; % choose from {'try', 'retry', 'load'}

if strcmp(mode, 'try')    
    [ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, done] = ...
        prep_world_coord_control_points(woco_image_file);
    
elseif strcmp(mode, 'retry')
    load(ctrl_pts_file);
    [ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, done] = ...
        prep_world_coord_control_points(...
            woco_image_file, ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, done);

elseif strcmp(mode, 'load')
    load(ctrl_pts_file);    

end

% save, take care to avoid accidental overwriting
if exist(ctrl_pts_file, 'file') == 2
    button = questdlg('Control points file exists! Overwrite it?', 'WARNING', 'Yes', 'No', 'No');
    if strcmp(button, 'Yes')
        save(ctrl_pts_file, 'ctrl_xw', 'ctrl_yw', 'ctrl_xp', 'ctrl_yp', 'done'); % precious!
    end
end


%% Rectify and crop 

% parameters
crop_xw = [-0.5, 0.65];
crop_yw = [ 0.002, 0.150]; % trim bottom to crop mylar out

% actions
[woco, xw, yw] = ...
    prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, ...
        imread(woco_image_file), true);

[rgb, ~, ~] = ...
    prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, ...
        imread(step_image_file), true);

%% Define a manual mask

mask_manual = prep_mask_manual(rgb);

%% Apply automatic masking model

% parameters
hue_lim = [0.02, 0.2];
value_lim = [0.15, 1.0];
entropy_lim = [0.6, 1];
entropy_len = 11;
morph_open_rad = 10;
morph_erode_rad = 2;

% actions
mask_auto = prep_mask_auto(rgb, hue_lim, value_lim, entropy_lim, entropy_len, ...
                    morph_open_rad, morph_erode_rad, true, true);

%% adaptive histogram equalization

% parameters
eql_len = 31;

% actions
eql = prep_intensity(rgb, mask_manual & mask_auto, eql_len, true, true);

%% run test analysis

prep_series(test_result_file, image_path, test_image_names, ctrl_xw, ...
    ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, hue_lim, ...
    value_lim, entropy_lim, entropy_len, morph_open_rad, ...
    morph_erode_rad, eql_len, xw, yw, mask_manual)

%% review results for test analysis

test_img = ncread(test_result_file, 'img');
test_img_rgb = ncread(test_result_file, 'img_rgb');
test_mask_auto = ncread(test_result_file, 'mask_auto');
test_mask_manual = ncread(test_result_file, 'mask_manual');
hf = figure;
for ii = 1:num_test_images
    hf.Name = sprintf('Test Analysis Results: %s', test_image_names{ii});
    subplot(3,1,1)
    imshow(test_img_rgb(:,:,:,ii));
    subplot(3,1,2)
    imshow(test_mask_manual & test_mask_auto(:,:,ii))
    subplot(3,1,3)
    imshow(test_img(:,:,ii));
    pause
end

%% Save parameters for batch processing

save(param_file, 'ctrl_xw', 'ctrl_yw', 'ctrl_xp', 'ctrl_yp', 'crop_xw', ...
    'crop_yw', 'hue_lim', 'value_lim', 'entropy_lim', 'entropy_len', ...
    'morph_open_rad', 'morph_erode_rad', 'eql_len', 'xw', 'yw', ...
    'mask_manual', 'result_file', 'image_path', 'image_names'); 
fprintf('Saved: %s\n', param_file);
