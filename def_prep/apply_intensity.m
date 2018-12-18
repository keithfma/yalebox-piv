function [] = apply_intensity(param_file, rgb, mask)
% function [] = apply_intensity(param_file, rgb, mask)
% 
% Apply histogram equalization using parameters specified in param_file
%
% Arguments:
%   param_file: string, path to parameter definition JSON file
%   img: 3D matrix, RGB color image OR 2D matrix, grayscale image
%   mask: 2D matrix, double, TRUE where there is sand, FALSE elsewhere
% %

param = loadjson(param_file, 'SimplifyCell', 1);

prep_intensity(...
    rgb, ...
    mask, ...
    param.intensity.eql_len.value, ...
    true, ...
    true);
