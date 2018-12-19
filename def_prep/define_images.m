function [] = define_images(param_file)
% function [] = define_images(param_file)
% 
% Define prep parameters in "images" section and write results to file
% 
% Arguments:
%   param_file: string, path to parameter definition JSON file
% %

param = loadjson(param_file, 'SimplifyCell', 1);

% path to images directory
path = uigetdir(...
    param.images.path.value, 'Select directory containing image files');
path = strip(path, 'right', filesep);

% world coordinate image
[woco_file, woco_path] = uigetfile('*.*', 'Select world coordinate image', ...
    fullfile(param.images.path.value, param.images.woco_file.value));
woco_path = strip(woco_path, 'right', filesep);
assert(strcmp(woco_path, path), 'Expect image to be in image directory');

% test image to use for parameter definition in this script
[test_file, test_path] = uigetfile('*.*', 'Select test image', ...
    fullfile(param.images.path.value, param.images.test_file.value));
test_path = strip(test_path, 'right', filesep);
assert(strcmp(test_path, path), 'Expect image to be in selected image directory');

% all experiment images
[exp_files, exp_path] = uigetfile(...
    '*.*', 'Select all experiment images', path, 'MultiSelect', 'on');
exp_path = strip(exp_path, 'right', filesep);
assert(strcmp(exp_path, path), 'Expect images to be in selected image directory');

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write "images" parameters to %s?', param_file);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.images.path.value = path;
    param.images.woco_file.value = woco_file;
    param.images.test_file.value = test_file;
    param.images.exp_files.value = exp_files;
    savejson('', param, 'Filename', param_file, 'SingletArray', 0);
    fprintf('"images" parameters written to file: %s\n', param_file);
else
    fprintf('"images" parameters NOT written to file\n');
end
