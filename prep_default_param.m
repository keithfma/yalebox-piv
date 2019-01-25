function pp = prep_default_param(filename)
% function pp = prep_default_param(filename)
%
% Save image pre-processing parameter file to specified path
% %

update_path('util');

[~, ~, ext] = fileparts(filename);
assert(strcmp(ext, '.json'), 'Filename must have extension .json');

% test ----------

pp.test.num_images.help = "int, number of images to include for test run";
pp.test.num_images.value = 20;

pp.test.test_file.help = "string, step image filename (without directory) to use for testing"; 
pp.test.test_file.value = "";

% images ---------

pp.images.path.help = "string, directory containing all image files";
pp.images.path.value = ""; 

pp.images.woco_file.help = "string, world coordinate image filename (without directory)";
pp.images.woco_file.value = "";

pp.images.exp_files.help = "list of strings, all experiment image filenames (without directory), assumed to be ordered by time";
pp.images.exp_files.value = [];

% woco ----------

pp.woco.xw.help = "1D vector, x-position of control points in world coordinates";
pp.woco.xw.value = [];

pp.woco.yw.help = "1D vector, x-position of control points in world coordinates";
pp.woco.yw.value = [];

pp.woco.xp.help = "1D vector, x-position of control points in pixel coordinates";
pp.woco.xp.value = [];

pp.woco.yp.help = "1D vector, x-position of control points in pixel coordinates";
pp.woco.yp.value = [];

% crop and rectify ----------

pp.crop.xlim.help = "2-element array, minimum and maximum x coordinates for crop in world coordinate units (m)";
pp.crop.xim.value = [-0.6, 0.6];

pp.crop.ylim.help = "2-element array, minimum and maximum y coordinates for crop in world coordinate units (m)";
pp.crop.ylim.value = [0.0, 0.2];

pp.crop.npts.help = "int, number of points to use in local weighted mean calculation for image warping";
pp.crop.npts.value = 10;

% manual masking ----------

pp.mask_manual.poly.help = "2D array, vertices of mask polygons, x-coords in row 1 and y-coords in row 2, polygons separated by NaN";
pp.mask_manual.poly.value = [[]];

% automatic masking ----------

pp.mask_auto.hue_lim.help = "2-element vector, double, range [0, 1]. [minimum, maximum] HSV 'hue' included as sand in the mask"; 
pp.mask_auto.hue_lim.value = [0.01, 0.30];

pp.mask_auto.value_lim.help = "2-element vector, double, range [0,1]. [minimum, maximum] HSV 'value' included as sand in the mask";
pp.mask_auto.value_lim.value = [0.07, 0.6];

pp.mask_auto.entropy_lim.help = "2-element vector, double, range [0, 1]. [minimum, maximum] entropy included as sand in the mask"; 
pp.mask_auto.entropy_lim.value = [0.55, 1.00]; 

pp.mask_auto.entropy_len.help = "scalar, integer, window size in pixels for entropy filter"; 
pp.mask_auto.entropy_len.value = 11;

pp.mask_auto.morph_open_rad.help = "scalar, double, radius of disk structuring element used in mophological opening filter";
pp.mask_auto.morph_open_rad.value =  10;

pp.mask_auto.morph_erode_rad.help = "scalar, double, radius of disk structuring element used in mophological erosion filter";
pp.mask_auto.morph_erode_rad.value = 3;

pp.intensity.eql_len.help = "scalar, integer, odd. Side length (in pixels) for the local neighborhood used to compute the transform for each pixel";
pp.intensity.eql_len.value = 41;

% save as file ---

save_param(pp, filename);
