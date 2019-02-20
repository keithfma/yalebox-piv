function pp = post_default_param(filename)
% function pp = post_default_param(filename)
%
% Save PIV post-processing parameter file to specified path
% %

update_path('util');

[~, ~, ext] = fileparts(filename);
assert(strcmp(ext, '.json'), 'Filename must have extension .json');

pp.notes.help = 'string, any notes to be included in the results file for posterity';
pp.notes.value = '';

pp.pro_bbox.help = '4-element position vector [xmin, ymin, width, height] for the bounding box used to estimate pro-side (mylar) displacement from PIV data results';
pp.pro_bbox.value = [NaN, NaN, NaN, NaN];

pp.retro_bbox.help = '4-element position vector [xmin, ymin, width, height] for the bounding box used to estimate retro-side (mylar) displacement from PIV data results';
pp.retro_bbox.value = [NaN, NaN, NaN, NaN];

pp.pad_method.help = 'string, selects which padding method to use prior to calculating spatial derivatives';
pp.pad_method.value = 'nearest';

save_param(pp, filename);
