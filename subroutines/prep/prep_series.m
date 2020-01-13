function [] = prep_series(result_file, image_path, image_names, image_view, ...
                  ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, ...
                  mask_poly, hue_lim, value_lim, entropy_lim, entropy_len, ...
                  equalize_len, pad_num_rows, pad_num_cols, fill_skin_min, ...
                  fill_skin_max, fill_bnd_smooth_window, fill_mirror, notes)
% function [] = prep_series(result_file, image_path, image_names, image_view, ...
%                   ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, ...
%                   mask_poly, hue_lim, value_lim, entropy_lim, entropy_len, ...
%                   equalize_len, pad_num_rows, pad_num_cols, fill_skin_min, ...
%                   fill_skin_max, fill_bnd_smooth_window, fill_mirror, notes)
% 
% Create PIV input file for a given image series. Reads in the images,
% rectifies and crops, masks, corrects illumination, and saves the results
% and metadata in a MAT file.
%
% Arguments:
% 
%   result_file = String, filename of the MAT file to be created. 
% 
%   image_path = String, path to folder containing the images.
%
%   image_names = Cell array of strings, cells must contain filenames for
%       successive images in the experiment image series. 
%
%   image_view = String, specify 'top' or 'side' view of experiment
%
%   ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp = Control points, as defined by
%       prep_world_coord_control_pts() for details
%
%   crop_xw, crop_yw = Image crop limits in world coordinates, see
%       prep_rectify_and_crop() for details
%   
%   mask_poly = 2D array, vertices of mask polygons, x-coords in row 1 and
%       y-coords in row 2, polygons separated by NaN
% 
%   hue_lim, value_lim, entropy_lim = 2-element vectors, threshold limits for
%       prep_mask_auto()
%
%   entropy_len = Size of entropy filter window, see prep_mask_auto()
%
%   equalize_len: Scalar, integer, odd. Side length (in pixels) for the local
%       neighborhood used to compute the transform for each pixel.
%
%   pad_num_rows, pad_num_cols: Scalar, integers, number of rows/cols to pad out on both sides
%
%   fill_skin_min, fill_skin_max: Scalar, integers, thickness of skin layer used in padding, see
%       prep_fill() for details.
%
%   fill_bnd_smooth_window: Integer, smoothing window to apply to boundary prior to fill, see
%       prep_fill() for details.
%
%   fill_mirror: bool, set true to mirror skin layer, or false to repeat it
%
%   notes: String, notes to be included in output MAT-file as a global
%       attribute. default = ''
% %

fprintf('%s: preprocessing %d images\n\n', mfilename, numel(image_names));

% load dependencies
update_path('prep', 'util');

% set defaults
narginchk(19, 20);
if nargin < 20; notes = ''; end

% check for sane arguments (pass-through arguments are checked in subroutines)
validateattributes(result_file, {'char'}, {'vector'});
validateattributes(image_path, {'char'}, {'vector'});
validateattributes(image_names, {'cell'}, {'vector'});
assert(any(strcmp(image_view, {'side', 'top'})), 'Invalid value for "image_view"');
[~, ~, result_file_ext] = fileparts(result_file);
assert(strcmp('.mat', result_file_ext), 'Output file must be .mat');

% check that all images exist and have the expected size and type
img_info = imfinfo(fullfile(image_path, image_names{1}));
img_nrow = img_info.Height;
img_ncol = img_info.Width;
for i = 1:numel(image_names)
    try
        this_file = fullfile(image_path, image_names{i});
        img_info = imfinfo(this_file);
        assert(img_info.Width == img_ncol && img_info.Height == img_nrow, ...
            sprintf('incorrect dimensions in image %s', this_file));
        assert(img_info.BitDepth == 24, ...
            sprintf('incorrect bit depth in image %s', this_file));
    catch
        error('unable to read image %i: %s', i, this_file);
    end
end

% create file, fail if exists
assert(exist(result_file, 'file') == 0, ...
    'Output file exists, either make space or choose another filename');
result = matfile(result_file, 'Writable', true);

% get coordinate vectors, some size parameters, and manual mask array
% note: created during prep steps, so prep a fake image)
fprintf('%s: get cordinate vectors and manual mask array\n', mfilename);
raw = zeros(img_nrow, img_ncol, 3, 'uint8'); 
[img, xw, yw] = prep_rectify_and_crop(...
    ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, raw);
[mask_manual, ~] = prep_mask_manual(img, mask_poly);
nx = numel(xw);
ny = numel(yw);
num_image = numel(image_names);

% define output metadata
meta = struct();
meta.notes = notes;
meta.version = get_version();
meta.view = image_view;

meta.args.ctrl_xw = ctrl_xw;
meta.args.ctrl_yw = ctrl_yw;
meta.args.ctrl_xp = ctrl_xp;
meta.args.ctrl_yp = ctrl_yp;
meta.args.crop_xw = crop_xw;
meta.args.crop_yw = crop_yw;
meta.args.mask_poly = mask_poly;
meta.args.hue_lim = hue_lim;
meta.args.value_lim = value_lim;
meta.args.entropy_lim = entropy_lim;
meta.args.entropy_len = entropy_len;
meta.args.equalize_len = equalize_len;
meta.args.pad_num_rows = pad_num_rows;
meta.args.pad_num_cols = pad_num_cols;
meta.args.fill_skin_min = fill_skin_min;
meta.args.fill_skin_max = fill_skin_max;
meta.args.fill_bnd_smooth_window = fill_bnd_smooth_window;
meta.args.fill_mirror = fill_mirror;

meta.x.name = 'x';
meta.x.long_name = 'horizontal position';
meta.x.dimensions = {};  % is coordinate axis
meta.x.notes = 'coordinate axis';
meta.x.units = 'meters';

meta.y.name = 'y';
meta.y.long_name = 'vertical position';
meta.y.notes = 'coordinate axis';
meta.y.dimensions = {};  % is coordinate axis
meta.y.units = 'meters';

meta.band.name = 'band';
meta.band.long_name = 'color band (RGB)';
meta.band.notes = 'coordinate axis';
meta.band.dimensions = {}; % is coordinate axis
meta.band.units = NaN;

meta.step.name = 'step';
meta.step.long_name = 'step number';
meta.step.notes = 'coordinate axis';
meta.step.dimensions = {};  % is coordinate axis
meta.step.units = '1';

meta.x_ext.name = 'x_ext';
meta.x_ext.long_name = 'horizonal position, extended';
meta.x_ext.notes = 'coordinate axis';
meta.x_ext.dimensions = {};
meta.x_ext.units = 'meters';

meta.y_ext.name = 'y_ext';
meta.y_ext.long_name = 'vertical position, extended';
meta.y_ext.notes = 'coordinate axis';
meta.y_ext.dimensions = {};
meta.y_ext.units = 'meters';

meta.img.name = 'img'; 
meta.img.long_name = 'rectified rgb image';
meta.img.notes = '';
meta.img.dimensions = {'y', 'x', 'band', 'step'};
meta.img.units = '24-bit color';

meta.img_eql.name = 'img_ext'; 
meta.img_eql.long_name = 'extended grayscale image';
meta.img_eql.notes = '';
meta.img_eql.dimensions = {'y_ext', 'x_ext', 'step'};
meta.img_eql.units = 'normalized';

meta.mask_auto.name = 'mask';
meta.mask_auto.long_name = 'sand mask';
meta.mask_auto.notes = '';
meta.mask_auto.dimensions = {'y', 'x', 'step'};
meta.mask_auto.units = 'boolean';

result.meta = meta;

% allocate output arrays
% note: store results in memory and write to file at once because
%   incremental writes to matfile slow to a crawl, allocating the variables
%   at the outset ensures that everything fits in memory
imgs = zeros(ny, nx, 3, num_image, 'uint8');
masks = false(ny, nx, num_image);

% keep track of the maximum extent of the mask
min_row = inf;  % set initial values so they are sure to be overridden
max_row = -inf;
min_col = inf;
max_col = -inf;

% read, rectify, crop, and mask all images
% note: will crop out empty space before pad, equalize and fill
parfor ii = 1:num_image
    
    this_file = fullfile(image_path, image_names{ii});
    img = imread(this_file);
     
    fprintf('\n%s: read, rectify, crop, and mask image %s\n', mfilename, this_file);
    
    img = prep_rectify_and_crop(...
        ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, img);

    mask_auto = prep_mask_auto(...
        img, hue_lim, value_lim, entropy_lim, entropy_len, image_view);
    mask = logical(mask_auto & mask_manual);
    
    [mask_rows, mask_cols] = ind2sub(size(mask), find(mask));
    min_row = min(min_row, min(mask_rows));
    max_row = max(max_row, max(mask_rows));
    min_col = min(min_col, min(mask_cols));
    max_col = max(max_col, max(mask_cols));    
    
    imgs(:, :, :, ii) = uint8(img);
    masks(:, :, ii) = mask;
    
end

% reduce size of data arrays
nx = max_col - min_col + 1;  % size of reduced data arrays
ny = max_row - min_row + 1;
xw = xw(min_col:max_col);
yw = yw(min_row:max_row);

imgs = imgs(min_row:max_row, min_col:max_col, :, :);
masks = masks(min_row:max_row, min_col:max_col, :);

% get coordinates and size parameters for extended arrays
[xw_ext, yw_ext, ~, ~] = prep_pad(xw, yw, ones(ny, nx), true(ny, nx), pad_num_rows, pad_num_cols);
nx_ext = numel(xw_ext);
ny_ext = numel(yw_ext);

% allocate extended arrays
imgs_ext = zeros(ny_ext, nx_ext, num_image, 'single');
masks_ext = false(ny_ext, nx_ext, num_image);

% equalize, pad, and fill all images
% note: will crop out empty space before pad, equalize and fill
parfor ii = 1:num_image
    
    img = imgs(:, :, :, ii);
    mask = masks(:, :, ii);
     
    fprintf('\n%s: equalize, pad, and fill image %i\n', mfilename, ii);
    
    img = prep_grayscale(img);
    
    [~, ~, img_ext, mask_ext] = prep_pad(...
        xw, yw, img, mask, pad_num_rows, pad_num_cols);
    
    [img_ext, mask_ext] = prep_fill(...
        img_ext, mask_ext, fill_skin_min, fill_skin_max, fill_bnd_smooth_window);
    
    img_ext = prep_equalize(img_ext, mask_ext, equalize_len);
    
    imgs_ext(:, :, ii) = single(img_ext);
    masks_ext(:, :, ii) = mask_ext;
    
end

% copy data from memory to results file
% note: here we trim image data to maximum extent of mask i.e., exclude
%   pixels that are always sand to reduce compute during PIV
fprintf('\n%s: copy from memory to results file: %s\n', mfilename, result_file);
fprintf('%s: images array size = %d x %d x 3 x %d\n', mfilename, ny, nx, num_image);
fprintf('%s: masks array size = %d x %d x %d\n', mfilename, ny, nx, num_image);
fprintf('%s: extended images array size = %d x %d x %d\n', mfilename, ny_ext, nx_ext, num_image);
fprintf('%s: extended masks array size = %d x %d x %d\n', mfilename, ny_ext, nx_ext, num_image);

result.x = xw;
result.y = yw;
result.x_ext = xw_ext;
result.y_ext = yw_ext;
result.band = 'RGB';
result.step = 0:(num_image - 1);

result.img = imgs;
result.mask = masks;

result.img_ext = imgs_ext;
result.mask_ext = masks_ext;
