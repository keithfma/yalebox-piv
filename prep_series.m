function [] = prep_series(result_file, image_path, image_names, image_view, ...
                  ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, fit_npts, ...
                  hue_lim, value_lim, entropy_lim, entropy_len, ...
                  morph_open_rad, morph_erode_rad, eql_len, xw, yw, mask_manual)
% function [] = prep_series(result_file, image_path, image_names, image_view, ...
%                   ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, fit_npts, ...
%                   hue_lim, value_lim, entropy_lim, entropy_len, ...
%                   morph_open_rad, morph_erode_rad, eql_len, xw, yw, mask_manual)
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
%       prep_world_coord_control_pts()
%
%   crop_xw, crop_yw = Image crop limits in world coordinates, see
%       prep_rectify_and_crop()
%
%   fit_npts = Local neighborhood for image rectification, see prep_rectify_and_crop()
%
%   hue_lim, value_lim, entropy_lim = 2-element vectors, threshold limits for
%       prep_mask_auto()
%
%   entropy_len = Size of entropy filter window, see prep_mask_auto()
%
%   morph_open_rad, morph_erode_rad = Scalar integers, structuring element
%       radius for morphological filters in prep_mask_auto()
%
%   eql_len = Neighborhood size for adaptive equalization, see prep_intensity()
%
%   xw, yw = Coordinate vectors for rectified/cropped image, in meters, see
%       prep_rectify_and_crop()
%
%   mask_manual = Output argument from prep_mask_manual()
% %

% load dependencies
update_path('prep', 'util');

% check for sane arguments (pass-through arguments are checked in subroutines)
assert(nargin == 21);
validateattributes(result_file, {'char'}, {'vector'});
validateattributes(image_path, {'char'}, {'vector'});
validateattributes(image_names, {'cell'}, {'vector'});
assert(any(strcmp(image_view, {'side', 'top'})), 'Invalid value for "image_view"');

[~, ~, result_file_ext] = fileparts(result_file);
assert(strcmp('.mat', result_file_ext), 'Output file must be .mat');

% get some size parameters
nx = numel(xw);
ny = numel(yw);
num_image = numel(image_names);

% check that all images exist and have the expected size and type
info = imfinfo(fullfile(image_path, image_names{1}));
raw_nrow = info.Height;
raw_ncol = info.Width;
for i = 1:num_image
    try
        this_file = fullfile(image_path, image_names{i});
        info = imfinfo(this_file);
        assert(info.Width == raw_ncol && info.Height == raw_nrow, ...
            sprintf('incorrect dimensions in image %s', this_file));
        assert(info.BitDepth == 24, ...
            sprintf('incorrect bit depth in image %s', this_file));
    catch
        error('unable to read image %i: %s', i, this_file);
    end
end

% initialize output file -------------------------------------------------

% create file, fail if exists
assert(exist(result_file, 'file') == 0, ...
    'Output file exists, either make space or choose another filename');
result = matfile(result_file, 'Writable', true);

% define metadata
meta = struct();

meta.version = get_version();

meta.view = image_view;

meta.input.ctrl_xw = ctrl_xw;
meta.input.ctrl_yw = ctrl_yw;
meta.input.ctrl_xp = ctrl_xp;
meta.input.ctrl_yp = ctrl_yp;
meta.input.crop_xw = crop_xw;
meta.input.crop_yw = crop_yw;
meta.input.fit_npts = fit_npts;
meta.input.hue_lim = hue_lim;
meta.input.value_lim = value_lim;
meta.input.entropy_lim = entropy_lim;
meta.input.entropy_len = entropy_len;
meta.input.morph_open_rad = morph_open_rad;
meta.input.morph_erode_rad = morph_erode_rad;
meta.input.eql_len = eql_len;

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

meta.img_rgb.name = 'img_rgb'; 
meta.img_rgb.long_name = 'rectified rgb image';
meta.img_rgb.notes = '';
meta.img_rgb.dimensions = {'y', 'x', 'rgb', 'step'};
meta.img_rgb.units = '24-bit color';

meta.img.name = 'img';
meta.img.long_name = 'rectified normalized grayscale image';
meta.img.notes = '';
meta.img.dimensions = {'y', 'x', 'step'};
meta.img.units = '1';

meta.mask_auto.name = 'mask_auto';
meta.mask_auto.long_name = 'sand mask, automatic';
meta.mask_auto.notes = '';
meta.mask_auto.dimensions = {'y', 'x', 'step'};
meta.mask_auto.units = 'boolean';

meta.mask_manual.name = 'mask_manual';
meta.mask_manual.long_name = 'sand mask, manual';
meta.mask_manual.notes = 'user-defined, constant in time';
meta.mask_manual.dimensions = {'y', 'x'};
meta.mask_manual.units = 'boolean';

result.meta = meta;

% write constant variables, allocate space for other variables
result.x = xw;

result.y = yw;

result.step = 0:(num_image - 1);

result.img_rgb = uint8.empty(0, 0, 0, 0);
result.img_rgb(ny, nx, 3, num_image) = uint8(0);

result.img = single.empty(0, 0, 0);
result.img(ny, nx, num_image) = single(0); 

result.mask_auto = logical.empty(0, 0, 0);
result.mask_auto(ny, nx, num_image) = false;

result.mask_manual = mask_manual;

% loop over all images
for ii = 1:num_image
    
    % read in original image
    this_file = fullfile(image_path, image_names{ii});
    raw = imread(this_file);
     
    % update user (always verbose)
    fprintf('\n%s: %s\n', mfilename, this_file);
    
    % rectify and crop
    raw = prep_rectify_and_crop(...
        ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, raw, fit_npts);
     
    % compute automatic mask
    mask_auto = prep_mask_auto(...
        raw, hue_lim, value_lim, entropy_lim, entropy_len, ...
        morph_open_rad, morph_erode_rad);
    
    % equalize intensity
    img = prep_intensity(raw, mask_manual & mask_auto, eql_len);
     
    % save results   
    result.mask_auto(:, :, ii) = logical(mask_auto);
    result.img(:, :, ii) = single(img);
    result.img_rgb(:, :, :, ii) = uint8(raw);
   
end
