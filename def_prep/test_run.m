function [] = test_run(param_file, num_images)
% Run prep for several images in the series and preview results from the
% output netCDF file
%
% Arguments:
%   param_file: string, path to parameter definition JSON file
%   num_images: number of images to include in the test run
% %

param = loadjson(param_file, 'SimplifyCell', 1);

% select evenly-spaced files from exp file list
idx = round(linspace(1, length(param.images.exp_files.value), num_images));
param.images.exp_files.value = param.images.exp_files.value(idx);

% save temporary param file, register for cleanup on function completion
temp_param_file = tempname;
clean_param = onCleanup(@() delete(temp_param_file));
savejson('', param, 'Filename', temp_param_file, 'SingletArray', 0);

% call prep series
temp_result_file = tempname;
clean_result = onCleanup(@() delete(temp_result_file));
% TODO: implement this wrapper function - make any changes needed to
%   preserve backward compatibility.
prep_series_from_file(temp_param_file, temp_result_file);

% plot results for each image in sequence
test_img = ncread(temp_result_file, 'img');
test_img_rgb = ncread(temp_result_file, 'img_rgb');
test_mask_auto = ncread(temp_result_file, 'mask_auto');
test_mask_manual = ncread(temp_result_file, 'mask_manual');

hf = figure;
for ii = 1:num_test_images
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