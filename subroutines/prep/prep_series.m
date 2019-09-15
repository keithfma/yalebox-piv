function [] = prep_series(result_file, image_path, image_names, image_view, ...
                  ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, ...
                  mask_poly, hue_lim, value_lim, entropy_lim, entropy_len, ...
                  equalize_len, notes)
% function [] = prep_series(result_file, image_path, image_names, image_view, ...
%                   ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, ...
%                   mask_poly, hue_lim, value_lim, entropy_lim, entropy_len, ...
%                   equalize_len, notes)
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
%   equalize_len_len: Scalar, integer, odd. Side length (in pixels) for the local
%       neighborhood used to compute the transform for each pixel.
%
%   notes: String, notes to be included in output MAT-file as a global
%       attribute. default = ''
% %

fprintf('%s: preprocessing %d images\n\n', mfilename, numel(image_names));

% load dependencies
update_path('prep', 'util');

% set defaults
narginchk(16, 17);
if nargin < 17; notes = ''; end

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

% create temporary output file, to hold results before crop to mask extent
% note: needed to avoid loading whole data arrays to memory during crop
% note: registers a function to cleanup on exit, no matter what
tmp_name = [tempname(), '.mat'];
tmp = matfile(tmp_name, 'Writable', true);
cleaner = onCleanup(@() delete(tmp_name));

% get coordinate vectors and manual mask array
% note: created during prep steps, so prep a fake image)
fprintf('%s: get cordinate vectors and manual mask array\n', mfilename);
raw = zeros(img_nrow, img_ncol, 3, 'uint8'); 
[img, xw, yw] = prep_rectify_and_crop(...
    ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, raw);
[mask_manual, ~] = prep_mask_manual(img, mask_poly);

% get some size parameters
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

meta.img.name = 'img'; 
meta.img.long_name = 'rectified rgb image';
meta.img.notes = '';
meta.img.dimensions = {'y', 'x', 'band', 'step'};
meta.img.units = '24-bit color';

meta.img_eql.name = 'img_eql'; 
meta.img_eql.long_name = 'histogram-equalized grayscale image';
meta.img_eql.notes = '';
meta.img_eql.dimensions = {'y', 'x', 'step'};
meta.img_eql.units = '1';

meta.mask_auto.name = 'mask';
meta.mask_auto.long_name = 'sand mask';
meta.mask_auto.notes = '';
meta.mask_auto.dimensions = {'y', 'x', 'step'};
meta.mask_auto.units = 'boolean';

tmp.meta = meta;

% write constant variables, allocate space for other variables
tmp.x = xw;
tmp.y = yw;
tmp.band = 'RGB';
tmp.step = 0:(num_image - 1);
allocate(tmp, 'img', 'uint8', [ny, nx, 3, num_image]);
allocate(tmp, 'img_eql', 'single', [ny, nx, num_image]);
allocate(tmp, 'mask', 'logical', [ny, nx, num_image]);

% keep track of the maximum extent of the mask
min_row = inf;  % set initial values so they are sure to be overridden
max_row = -inf;
min_col = inf;
max_col = -inf;

% loop over all images
for ii = 1:num_image
    
    this_file = fullfile(image_path, image_names{ii});
    img = imread(this_file);
     
    fprintf('\n%s: preprocess image %s\n', mfilename, this_file);
    
    img = prep_rectify_and_crop(...
        ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, img);

    mask_auto = prep_mask_auto(...
        img, hue_lim, value_lim, entropy_lim, entropy_len, image_view);
    mask = logical(mask_auto & mask_manual);
    
    img_eql = prep_equalize(img, mask, equalize_len);
    
    [mask_rows, mask_cols] = ind2sub(size(mask), find(mask));
    min_row = min(min_row, min(mask_rows));
    max_row = max(max_row, max(mask_rows));
    min_col = min(min_col, min(mask_cols));
    max_col = max(max_col, max(mask_cols));    
    
    tmp.img(:, :, :, ii) = uint8(img);
    tmp.img_eql(:, :, ii) = single(img_eql);
    tmp.mask(:, :, ii) = mask;

end
nx = max_col - min_col + 1;  % size of reduced data arrays
ny = max_row - min_row + 1;

% copy data from temporary file to results file
% note: here we trim image data to maximum extent of mask i.e., exclude
%   pixels that are always sand to reduce compute during PIV
% note: copy data in chunks to avoid blowing up memory
fprintf('\n%s: copy from temporary to final file: %s\n', mfilename, result_file);
fprintf('%s: images array size = %d x %d x 3 x %d\n', mfilename, ny, nx, num_image);
fprintf('%s: masks array size = %d x %d x %d\n', mfilename, ny, nx, num_image);

% compute chunk size that fits in constant memory limits
max_memory_bytes = 8*1024^3;  % 8 GB
image_bytes = ny*nx*3*1; % uint8s are 1 byte
image_eql_bytes = ny*nx*32;  % singles are 32 bytes
mask_bytes = ny*nx*1;  % logicals are 1 byte
total_bytes = image_bytes + image_eql_bytes + mask_bytes;
chunk_size = floor(max_memory_bytes/total_bytes);
fprintf('%s: copying in chunks of %d images\n', mfilename, chunk_size);

result.meta = tmp.meta;
result.x = tmp.x(1, min_col:max_col);
result.y = tmp.y(1, min_row:max_row);
result.band = tmp.band;
result.step = tmp.step;

allocate(result, 'img', 'uint8', [ny, nx, 3, num_image]);
allocate(result, 'img_eql', 'single', [ny, nx, num_image]);
allocate(result, 'mask', 'logical', [ny, nx, num_image]);

for min_img = 1:chunk_size:num_image
    max_img = min(min_img + chunk_size - 1, num_image);
    fprintf('%s: copy images %d-%d\n', mfilename, min_img, max_img);
    result.img(:, :, :, min_img:max_img) = ...
        tmp.img(min_row:max_row, min_col:max_col, :, min_img:max_img);
    result.img_eql(:, :, min_img:max_img) = ...
        tmp.img_eql(min_row:max_row, min_col:max_col, min_img:max_img);
    result.mask(:, :, min_img:max_img) = ...
        tmp.mask(min_row:max_row, min_col:max_col, min_img:max_img);
end
