function [] = prep_series(result_file, image_path, image_names, image_view, ...
                  ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, ...
                  mask_poly, hue_lim, value_lim, entropy_lim, entropy_len, ...
                  morph_open_rad, morph_erode_rad, notes)
% function [] = prep_series(result_file, image_path, image_names, image_view, ...
%                   ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, ...
%                   mask_poly, hue_lim, value_lim, entropy_lim, entropy_len...
%                   morph_open_rad, morph_erode_rad, notes)
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
%   morph_open_rad, morph_erode_rad = Scalar integers, structuring element
%       radius for morphological filters in prep_mask_auto()
%
%   notes: String, notes to be included in output MAT-file as a global
%       attribute. default = ''
% %

% load dependencies
update_path('prep', 'util');

% set defaults
narginchk(17, 18);
if nargin < 18; notes = ''; end

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

% get coordinate vectors and manual mask array
% note: created during prep steps, so prep a fake image)
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
meta.args.morph_open_rad = morph_open_rad;
meta.args.morph_erode_rad = morph_erode_rad;

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

meta.step.name = 'step';
meta.step.long_name = 'step number';
meta.step.notes = 'coordinate axis';
meta.step.dimensions = {};  % is coordinate axis
meta.step.units = '1';

meta.img.name = 'img'; 
meta.img.long_name = 'rectified rgb image';
meta.img.notes = '';
meta.img.dimensions = {'y', 'x', 'rgb', 'step'};
meta.img.units = '24-bit color';

meta.mask_auto.name = 'mask';
meta.mask_auto.long_name = 'sand mask';
meta.mask_auto.notes = '';
meta.mask_auto.dimensions = {'y', 'x', 'step'};
meta.mask_auto.units = 'boolean';

result.meta = meta;

% write constant variables, allocate space for other variables
result.x = xw;
result.y = yw;
result.step = 0:(num_image - 1);
allocate(result, 'img', 'uint8', [ny, nx, 3, num_image]);
allocate(result, 'mask', 'logical', [ny, nx, num_image]);

% loop over all images
for ii = 1:num_image
    
    this_file = fullfile(image_path, image_names{ii});
    img = imread(this_file);
     
    fprintf('\n%s: %s\n', mfilename, this_file);
    
    img = prep_rectify_and_crop(...
        ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, img);

    mask_auto = prep_mask_auto(...
        img, hue_lim, value_lim, entropy_lim, entropy_len, ...
        morph_open_rad, morph_erode_rad);
      
    result.mask(:, :, ii) = logical(mask_auto & mask_manual);
    result.img(:, :, :, ii) = uint8(img);

end
