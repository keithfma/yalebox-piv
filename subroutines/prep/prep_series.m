function [] = prep_series(result_file, image_path, image_names, image_view, ...
                  ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, ...
                  segment_scale, segment_sigma, segment_min_size, ...
                  train_features, train_labels, eql_len, notes)
% function [] = prep_series(result_file, image_path, image_names, image_view, ...
%                   ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, ...
%                   segment_scale, segment_sigma, segment_min_size, ...
%                   train_features, train_labels, eql_len, notes)
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
%   segment_scale, segment_sigma, segment_min_size: Image segmentation
%       parameters, see prep_mask_segment() for details
%   
%   train_features, train_labels: mask model training data, see
%       prep_mask_train() for details
% 
%   eql_len: histogram equalization kernel size, see prep_intensity() for
%       details
%
%   notes: String, notes to be included in output MAT-file as a global
%       attribute. default = ''
% %

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

% get coordinate vectors (created during prep steps, so prep a fake image)
img_fake = zeros(img_nrow, img_ncol, 3, 'uint8'); 
[~, xw, yw] = prep_rectify_and_crop(...
    ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, img_fake);

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
meta.args.segment_scale = segment_scale;
meta.args.segment_sigma = segment_sigma;
meta.args.segment_min_size = segment_min_size;
meta.args.train_features = train_features;
meta.args.train_labels = train_labels;

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

meta.mask.name = 'mask';
meta.mask.long_name = 'sand mask';
meta.mask.notes = 'true where pixel is sand, false elsewhere';
meta.mask.dimensions = {'y', 'x', 'step'};
meta.mask.units = 'boolean';

result.meta = meta;

% write constant variables, allocate space for other variables
result.x = xw;
result.y = yw;
result.step = 0:(num_image - 1);
allocate(result, 'img_rgb', 'uint8', [ny, nx, 3, num_image]);
allocate(result, 'img', 'single', [ny, nx, num_image]); 
allocate(result, 'mask', 'logical', [ny, nx, num_image]);

% train mask classifier model
% note: not stored to output avoid version hell
fprintf('%s: train mask classifier model\n', mfilename);
mask_classifier = prep_mask_train(train_features, train_labels);

% loop over all images
for ii = 1:num_image
    
    this_file = fullfile(image_path, image_names{ii});
    img_rgb = imread(this_file);
     
    fprintf('\n%s: %s\n', mfilename, this_file);
    
    img_rgb = prep_rectify_and_crop(...
        ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, img_rgb);

    img_segments = prep_mask_segments(...
        img_rgb, segment_scale, segment_sigma, segment_min_size);
    
    img_features = prep_mask_features(img_rgb, img_segments);
    
    img_mask = prep_mask_apply...
        (mask_classifier, img_features, img_segments);
    
    % note: this step is legacy, remove after testing multiband correlation
    img = prep_intensity(img_rgb, img_mask, eql_len);
     
    result.mask(:, :, ii) = logical(img_mask);
    result.img(:, :, ii) = single(img);
    result.img_rgb(:, :, :, ii) = uint8(img_rgb);
   
end
