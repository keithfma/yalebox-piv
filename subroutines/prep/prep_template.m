function pp = prep_template(filename)
% function pp = prep_template(filename)
%
% Generate and optionally save template image pre-processing parameters
% 
% Arguments:
% 
%   filename: optional string, path to output MAT file, if not provided,
%       parameters will be created but not saved
% 
% Returns:
%   
%   template parameter struct (always), and saves parameter struct to file
%       if filename is providede
% %

update_path('util');

if nargin > 0
    [~, ~, ext] = fileparts(filename);
    assert(strcmp(ext, '.mat'), 'Filename must have extension .mat');
end

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

% equalization -----------

pp.equalize_len.help =  'scalar, integer, odd, side length (in pixels) for the local neighborhood used to compute the transform for each pixel';
pp.equalize_len.value = 21;

% pad and fill -----------

pp.pad_num_rows.help = 'Number of rows of padding to add to top and bottom of all images, use to ensure extra displacement observations are available for gradient calculation';
pp.pad_num_rows.value = 100;

pp.pad_num_cols.help = 'Number of columns of padding to add to left and right of all images, use to ensure extra displacement observations are available for gradient calculation';
pp.pad_num_cols.value = 100;

pp.fill_skin_min.help = 'Integer, minimum thickness of skin layer to use when filling';
pp.fill_skin_min.value = 3;

pp.fill_skin_max.help = 'Integer, maximum thickness of skin layer to use when filling';
pp.fill_skin_max.value = 10;

pp.fill_bnd_smooth_window.help = 'width of the (loess) smoothing window to apply to the upper boundary before computing the skin layer, set to 0 to disable smoothing';
pp.fill_bnd_smooth_window.value = 10;

pp.fill_mirror.help = 'logical, set True to mirror the skin layer when filling, or false to repeat it';
pp.fill_mirror.value = false;

% save as file ---

if nargin > 0
    save_param(pp, filename);
end
