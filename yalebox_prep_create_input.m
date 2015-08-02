function [] = yalebox_prep_create_input(input_file, image_path, image_names, ...
                    x, y, x_scale, y_scale, x_offset, y_offset, mask_manual, ...
                    hue_lim, value_lim, entropy_lim, median_window, ...
                    entropy_window, opening_radius, histeq_width)
% 
% Create PIV input file for a given image series. Reads in the images,
% performs masking and color correction, and saves the results and metadata
% in a netCDF file.
%
% Arguments:
% 
% input_file = String, filename of the netCDF input file to be created. 
%
% image_path = String, path to folder containing the images.
%
% images_names = Cell array of strings, cells must contain filenames for
%   successive images in the experiment image series. 
%
% x, y, x_scale, y_scale, x_offset, y_offset = Output arguments from
%   yalebox_prep_world_coord().
%
% mask_manual = Output argument from yalebox_prep_mask_manual()
%
% hue_lim, value_lim, entropy_lim, median_window, entropy_window,
%   opening_radius = Select input arguments from yalebox_prep_mask_auto()
% 
% histeq_width = Select input argument from yalebox_prep_intensity()
%
% PIV input netCDF format:
%   groups: root, preprocess, input
%   dimensions: root/x, root/y, root/step
%   variables: preprocess/input_file, mask_auto, mask_manual, input/intensity
%   attributes: preprocess/* for all preprocessing parameters
%
% Keith Ma, July 2015

% check for sane arguments - only arguments that are unique to this
%   function are checked, those that are passed to the yalebox_prep_*
%   functions are checked in this procedues.
assert(isa(input_file, 'char'), ...
    'input_file is not a string.');
assert(isa(image_path, 'char') && exist(image_path, 'dir') == 7, ...
    'image_path is not a directory');
assert(isa(image_names, 'cell'), ...
    'image_names is not a cell array');

% check that all images exist and have the expected size and type
nimage = numel(image_names);
image_w = size(mask_manual, 2);
image_h = size(mask_manual, 1);
for i = 1:nimage
    try
        this_file = [image_path filesep image_names{i}];
        info = imfinfo(this_file);
        assert(info.Width == image_w && info.Height == image_h, ...
            sprintf('incorrect dimensions in image %s', this_file));
        assert(info.BitDepth == 24, ...
            sprintf('incorrect bit depth in image %s', this_file));
    catch 
        error(sprintf('Unable to read image %i: %s', i, this_file));
    end
end

% % check parameter struct members 
% 
% 
%        
% 
% 
% % check that all of the image files exist
% 
% % check that image and mask dimensions match the coordinate vectors
% 
% 
% % for i = 1:numel(image_files)
% %     assert(exist(image_files{i}, 'file') == 2, ...
% %         sprintf('image_files{%i} does not exist', i));
% %    
% %     info = iminfo(image_files{i});
% %     if i>1
% %         assert(info.Width == info_prev.Width && info.Height == info_prev.Height, ...
% %             sprintf('image_files{%i} is not the same size as the previous image', i));
% %     end
% %     info_prev = info;
% % end
% %         
%         
% % CHECK THE CONTENTS OF THE PARAM MATRIX
% 
% % init netcdf file
% % ADD PARAMETERS
% % ADD DIMENSIONS AND VARIABLES
% 
% % loop over all images
% 
% % read in original image
