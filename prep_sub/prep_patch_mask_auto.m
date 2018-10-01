function prep_patch_mask_auto(image_file, test_step, hue_lim, value_lim, ...
            entropy_lim, entropy_len, morph_open_rad, morph_erode_rad)
% function prep_patch_mask_auto(image_file, test_step, hue_lim, value_lim, ...
%             entropy_lim, entropy_len, morph_open_rad, morph_erode_rad)
% 
% Helper function to patch (i.e. overwrite) the manual_auto variable in an image
% series netCDF file (as produced by prep_series)
%
% NOTE: The field is *overwritten* in the input file, but you will be prompted
%   to confirm it first given the mask results for the test step
%
% Arguments:
%   image_file: String, image series netCDF as produced by prep_series
%   test_step: Scalar integer, index of step to display for mask definition
%   hue_lim, value_lim, etc: Automatic masking routine parameters, see
%       prep_mask_auto help documentation for details
%
% % Keith Ma

% check inputs
assert(exist(image_file, 'file') == 2);
validateattributes(test_step, {'numeric'}, {'scalar', 'positive', 'integer'});

% create and examine automatic mask for test step
img = ncread(image_file, 'img_rgb', [1, 1, 1, test_step], [inf, inf, 3, 1]);
prep_mask_auto(img, hue_lim, value_lim, entropy_lim, entropy_len, ...
    morph_open_rad, morph_erode_rad, true, true);

% confirm user intent to overwrite, then overwrite
button = questdlg(...
    sprintf('Do you *really* want to overwrite the "mask_auto" variable in %s?', image_file), ...
    mfilename, 'Yes', 'No', 'No');

if strcmp(button, 'Yes')
    
    % update attributes
    ncwriteatt(image_file, '/', 'prep_mask_auto hue_lim', hue_lim);
    ncwriteatt(image_file, '/', 'prep_mask_auto value_lim', value_lim);
    ncwriteatt(image_file, '/', 'prep_mask_auto entropy_lim', entropy_lim);
    ncwriteatt(image_file, '/', 'prep_mask_auto entropy_len', entropy_len);
    ncwriteatt(image_file, '/', 'prep_mask_auto morph_open_rad', morph_open_rad);
    ncwriteatt(image_file, '/', 'prep_mask_auto morph_erode_rad', morph_erode_rad);

    % update frames
    num_step = length(ncread(image_file, 'step'));
    for ii = 1:num_step
        fprintf('%d of %d\n', ii, num_step);
        img = ncread(image_file, 'img_rgb', [1, 1, 1, ii], [inf, inf, 3, 1]);
        mask = prep_mask_auto(img, hue_lim, value_lim, entropy_lim, ...
            entropy_len, morph_open_rad, morph_erode_rad, false, true);
        ncwrite(image_file, 'mask_auto', uint8(mask), [1, 1, ii]);
    end
       
    fprintf('%s: Completed\n', mfilename);

else
    fprintf('%s: Canceled\n', mfilename);
end
