function prep_patch_mask_manual(image_file, test_step)
% function prep_patch_mask_manual(original_file, patched_file, test_step)
% 
% Helper function to patch (i.e. overwrite) the manual_mask variable in an image
% series netCDF file (as produced by prep_series)
%
% NOTE: The field is *overwritten* in the input file, but you will be prompted
%   to confirm it first
%
% Arguments:
%   image_file: String, image series netCDF as produced by prep_series
%   test_step: Scalar integer, index of step to display for mask definition
%
% % Keith Ma

% check inputs
assert(exist(image_file, 'file') == 2);
validateattributes(test_step, {'numeric'}, {'scalar', 'positive', 'integer'});

% define mask, manually
img = ncread(image_file, 'img_rgb', [1, 1, 1, test_step], [inf, inf, 3, 1]);
mask = prep_mask_manual(img);

% confirm user intent to overwrite, then overwrite
button = questdlg(...
    sprintf('Do you *really* want to overwrite the "mask_manual" variable in %s?', image_file), ...
    mfilename, 'Yes', 'No', 'No');
if strcmp(button, 'Yes')
    fprintf('%s: Completed\n', mfilename);
    ncwrite(image_file, 'mask_manual', uint8(mask));
else
    fprintf('%s: Canceled\n', mfilename);
end
