function rgb = apply_crop(param_file)
% function rgb = apply_crop(param_file)
% 
% Apply rectification using control points, and crop specified by
% parameters, then display the results
% 
% Arguments:
%   param_file: string, path to parameter definition JSON file
% 
% Returns:
%   rgb: N x M x 3 matrix, cropped and rectified RGB image
% % 

param = loadjson(param_file, 'SimplifyCell', 1);

fprintf('Displaying rectifed & cropped woco image using current parameters\n');
[~, ~, ~] = ...1111
    prep_rectify_and_crop(...
        param.woco.xp.value, ...
        param.woco.yp.value, ...
        param.woco.xw.value, ...
        param.woco.yw.value, ...
        param.crop.xlim.value, ...
        param.crop.ylim.value, ...
        imread(fullfile(param.images.path.value, param.images.woco_file.value)), ...
        [], ...
        true);

fprintf('Displaying rectifed & cropped test image using current parameters\n');
[rgb, ~, ~] = ...
    prep_rectify_and_crop(...
        param.woco.xp.value, ...
        param.woco.yp.value, ...
        param.woco.xw.value, ...
        param.woco.yw.value, ...
        param.crop.xlim.value, ...
        param.crop.ylim.value, ...
        imread(fullfile(param.images.path.value, param.images.test_file.value)), ...,
        [], ...
        true);
