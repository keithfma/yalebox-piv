function pp = prep_template(filename)
% function pp = prep_template(filename)
%
% Save template for image pre-processing parameter file to specified path
% 
% Arguments:
% 
%   filename: string, path to output JSON file 
% %

update_path('util');

[~, ~, ext] = fileparts(filename);
assert(strcmp(ext, '.json'), 'Filename must have extension .json');

% general --------

pp.notes.help = 'string, any notes to be included in the results file for posterity';
pp.notes.value = '';

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

% masking ----------

pp.mask_train_sand.help = "Cell array of 2D arrays, vertices of ROI polygons for 'sand' class training data, each cell contains vertices for one polygon with x-coord in column 1 and y-coord in column 2";
pp.mask_train_sand.value = {};

pp.mask_train_other.help = "Cell array of 2D arrays, vertices of ROI polygons for 'other' class training data, each cell contains vertices for one polygon with x-coord in column 1 and y-coord in column 2";
pp.mask_train_other.value = {};

pp.mask_train_file.help = "string, filename (without directory) of image to use for mask model training data";
pp.mask_train_file.value = '';

pp.mask_model_type.help = "string, select mask model from a few options: 'tree' uses a simple decision tree (relatively fast for development), and 'forest' uses a slower, better random forest model.";
pp.mask_model_type.value = 'tree';

% obsolete?
pp.mask_poly.help = "2D array, vertices of mask polygons, x-coords in row 1 and y-coords in row 2, polygons separated by NaN";
pp.mask_poly.value = [[]];

% obsolete?
pp.mask_hue_lim.help = "2-element vector, double, range [0, 1]. [minimum, maximum] HSV 'hue' included as sand in the mask"; 
pp.mask_hue_lim.value = [0.01, 0.30];

% obsolete?
pp.mask_value_lim.help = "2-element vector, double, range [0,1]. [minimum, maximum] HSV 'value' included as sand in the mask";
pp.mask_value_lim.value = [0.07, 0.6];

% obsolete?
pp.mask_entropy_lim.help = "2-element vector, double, range [0, 1]. [minimum, maximum] entropy included as sand in the mask"; 
pp.mask_entropy_lim.value = [0.55, 1.00]; 

% obsolete?
pp.mask_entropy_len.help = "scalar, integer, window size in pixels for entropy filter"; 
pp.mask_entropy_len.value = 11;

% obsolete?
pp.mask_morph_open_rad.help = "scalar, double, radius of disk structuring element used in mophological opening filter";
pp.mask_morph_open_rad.value =  10;

% obsolete?
pp.mask_morph_erode_rad.help = "scalar, double, radius of disk structuring element used in mophological erosion filter";
pp.mask_morph_erode_rad.value = 3;

% color adjustment ----------

pp.intensity_eql_len.help = "scalar, integer, odd. Side length (in pixels) for the local neighborhood used to compute the transform for each pixel";
pp.intensity_eql_len.value = 41;

% save as file ---

save_param(pp, filename);
