function [] = define_images(param_file)
% function [] = define_images(param_file)
% 
% Define prep parameters in "images" section and write results to file
% 
% Arguments:
%   param_file: string, path to parameter definition JSON file
% %

param = loadjson(param_file, 'SimplifyCell', 1);

% world coordinate image
[woco_file_name, woco_file_path] = uigetfile(...
    '*.*', 'Select World Coordinate Image', param.images.woco_file.value);
woco_file = fullfile(woco_file_path, woco_file_name);

% test image to use for parameter definition in this script
[test_file_name, test_file_path] = uigetfile(...
    '*.*', 'Select Test Image', param.images.test_file.value);
test_file = fullfile(test_file_path, test_file_name);

% all experiment images
ini_exp_file = '';
if ~isempty(param.images.exp_files.value)
    ini_exp_file = param.images.exp_files.value{1};
end
[exp_files, exp_files_path] = uigetfile(...
    '*.*', 'Select All Experiment Images', ini_exp_file, 'MultiSelect', 'on');
for ii = 1:length(exp_files)
    exp_files{ii} = fullfile(exp_files_path, exp_files{ii});
end

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write "images" parameters to %s?', param_file);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.images.woco_file.value = woco_file;
    param.images.test_file.value = test_file;
    param.images.exp_files.value = exp_files;
    savejson('', param, 'Filename', param_file, 'SingletArray', 0);
    fprintf('"images" parameters written to file: %s\n', param_file);
else
    fprintf('"images" parameters NOT written to file\n');
end
