function [] = define_woco(param_file)
% function [] = define_woco(param_file)
% 
% Define parameters in "woco" section by creating or editing world
% coordinate control points interactively
% 
% Arguments:
%   param_file: string, path to parameter definition JSON file
% % 

param = loadjson(param_file, 'SimplifyCell', 1);

% sanity check
lengths = [length(param.woco.xw.value), ...
    length(param.woco.yw.value), ...
    length(param.woco.xp.value), ...
    length(param.woco.yp.value)];
if min(lengths) ~= max(lengths)
    error('Length of control point vectors do not all match');
end

% define control points, either from scratch or starting from existing
if max(lengths) > 0
    [xw, yw, xp, yp, ~] = prep_world_coord_control_points(...
        param.images.woco_file.value, ...
        param.woco.xw.value, ...
        param.woco.yw.value, ...
        param.woco.xp.value, ...
        param.woco.yp.value, ...
        true(size(param.woco.xw.value)));
else
    [xw, yw, xp, yp, ~] = prep_world_coord_control_points(param.images.woco_file.value);
end

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write "woco" parameters to %s?', param_file);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.woco.xw.value = xw;
    param.woco.yw.value = yw;
    param.woco.xp.value = xp;
    param.woco.yp.value = yp;
    savejson('', param, 'Filename', param_file, 'SingletArray', 0);
    fprintf('"woco" parameters written to file: %s\n', param_file);
else
    fprintf('"woco" parameters NOT written to file\n');
end
