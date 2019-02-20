function pp = prep_default_param(filename)
% function pp = prep_default_param(filename)
%
% Save image pre-processing parameter file to specified path
% %

update_path('util');

[~, ~, ext] = fileparts(filename);
assert(strcmp(ext, '.json'), 'Filename must have extension .json');

% images ---------

pp.view.help = "string, specify 'side' or 'top' image view";
pp.view.value = 'side';

pp.image_dir.help = "string, directory containing all image files";
pp.image_dir.value = ""; 

pp.woco_file.help = "string, world coordinate image filename (without directory)";
pp.woco_file.value = "";

pp.exp_files.help = "list of strings, all experiment image filenames (without directory), assumed to be ordered by time";
pp.exp_files.value = [];

% woco ----------

pp.woco_xw.help = "1D vector, x-position of control points in world coordinates";
pp.woco_xw.value = [];

pp.woco_yw.help = "1D vector, x-position of control points in world coordinates";
pp.woco_yw.value = [];

pp.woco_xp.help = "1D vector, x-position of control points in pixel coordinates";
pp.woco_xp.value = [];

pp.woco_yp.help = "1D vector, x-position of control points in pixel coordinates";
pp.woco_yp.value = [];

% crop and rectify ----------

pp.crop_xlim.help = "2-element array, minimum and maximum x coordinates for crop in world coordinate units (m)";
pp.crop_xlim.value = [-0.6, 0.6];

pp.crop_ylim.help = "2-element array, minimum and maximum y coordinates for crop in world coordinate units (m)";
pp.crop_ylim.value = [0.0, 0.2];

pp.crop_npts.help = "int, number of points to use in local weighted mean calculation for image warping";
pp.crop_npts.value = 10;

% masking ----------

pp.mask_poly.help = "2D array, vertices of mask polygons, x-coords in row 1 and y-coords in row 2, polygons separated by NaN";
pp.mask_poly.value = [[]];

pp.mask_hue_lim.help = "2-element vector, double, range [0, 1]. [minimum, maximum] HSV 'hue' included as sand in the mask"; 
pp.mask_hue_lim.value = [0.01, 0.30];

pp.mask_value_lim.help = "2-element vector, double, range [0,1]. [minimum, maximum] HSV 'value' included as sand in the mask";
pp.mask_value_lim.value = [0.07, 0.6];

pp.mask_entropy_lim.help = "2-element vector, double, range [0, 1]. [minimum, maximum] entropy included as sand in the mask"; 
pp.mask_entropy_lim.value = [0.55, 1.00]; 

pp.mask_entropy_len.help = "scalar, integer, window size in pixels for entropy filter"; 
pp.mask_entropy_len.value = 11;

pp.mask_morph_open_rad.help = "scalar, double, radius of disk structuring element used in mophological opening filter";
pp.mask_morph_open_rad.value =  10;

pp.mask_morph_erode_rad.help = "scalar, double, radius of disk structuring element used in mophological erosion filter";
pp.mask_morph_erode_rad.value = 3;

% color adjustment ----------

pp.intensity_eql_len.help = "scalar, integer, odd. Side length (in pixels) for the local neighborhood used to compute the transform for each pixel";
pp.intensity_eql_len.value = 41;

% save as file ---

save_param(pp, filename);