function [mrgb, mask] = define_mask_manual(param_file, rgb)
% function [mrgb, mask] = define_mask_manual(param_file, rgb)
% 
% Define parameters in "mask_manual" section by interactively editing
% 
% Arguments:
%   rgb: N x M x 3 matrix, cropped and rectified RGB image
%   param_file: string, path to parameter definition JSON file
% 
% Returns:
%   mrgb: N x M x 3 matrix, Rmrgb, GB image with manual mask applied to it
%   mask: N x M logical matrix
% % 

param = loadjson(param_file, 'SimplifyCell', 1);

% interactive mask creation
[mask, poly] = prep_mask_manual(...
    rgb, ...
    param.mask_manual.poly.value, ...
    true);

% create masked image for next step
mrgb = rgb;
mrgb(repmat(~mask, 1, 1, 3)) = 0;

% save results, take care to avoid accidental overwriting
prompt = sprintf('Write "mask_manual" parameters to %s?', param_file);
button = questdlg(prompt, 'WARNING', 'Yes', 'No', 'Yes');
if strcmp(button, 'Yes')
    param.mask_manual.poly.value = poly;
    savejson('', param, 'Filename', param_file, 'SingletArray', 0);
    fprintf('"mask_manual" parameters written to file: %s\n', param_file);
else
    fprintf('"mask_manual" parameters NOT written to file\n');
end
