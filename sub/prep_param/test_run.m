function [] = test_run(param_file)
% Run prep for several images in the series and preview results from the
% output netCDF file
%
% Arguments:
%   param_file: string, path to parameter definition JSON file
% %

param = loadjson(param_file, 'SimplifyCell', 1);

% select evenly-spaced files from exp file list
idx = round(linspace(...
    1, ...
    length(param.images.exp_files.value), ...
    param.test_run.num_images.value));
param.images.exp_files.value = param.images.exp_files.value(idx);

% save temporary param file, register for cleanup on function completion
temp_param_file = get_temp_file('json');
clean_param = onCleanup(@() delete(temp_param_file));
savejson('', param, 'Filename', temp_param_file, 'SingletArray', 0);

% call prep series
result_file = get_temp_file('nc');
fprintf('Test run results: %s\n', result_file);
prep_series_from_file(result_file, temp_param_file);

% plot results for each image in sequence
test_img = ncread(result_file, 'img');
test_img_rgb = ncread(result_file, 'img_rgb');
test_mask_auto = ncread(result_file, 'mask_auto');
test_mask_manual = ncread(result_file, 'mask_manual');

hf = figure;
for ii = 1:param.test_run.num_images.value
    hf.Name = sprintf('Test Analysis Results: %s', ...
        param.images.exp_files.value{ii});
    subplot(3,1,1)
    imshow(test_img_rgb(:,:,:,ii));
    subplot(3,1,2)
    imshow(test_mask_manual & test_mask_auto(:,:,ii))
    subplot(3,1,3)
    imshow(test_img(:,:,ii));
    pause
end