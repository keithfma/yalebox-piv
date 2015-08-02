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

% % check for sane arguments
% assert(isa(input_file, 'char'), ...
%     'input_file is not a string.');
% assert(isa(image_files, 'cell'), ...
%     'image_files is not a cell array');
% assert(isa(mask_manual, 'logical') && ismatrix(mask_manual), ...
%     'mask_manual is not a logical matrix');
% assert(isstruct(mask_auto),...
%     'mask_auto is not a struct .');
% assert(ismember('hue_lim', fieldnames(mask_auto)) && ...
%     ismember('value_lim', fieldnames(mask_auto)) && ...
%     ismember('entropy_lim', fieldnames(mask_auto)) && ...
%     ismember('median_window', fieldnames(mask_auto)) && ...
%     ismember('entropy_window', fieldnames(mask_auto)) && ...
%     ismember('opening_radius', fieldnames(mask_auto)), ...
%     'mask_auto does not contain all of the required parameters');
% assert(isstruct(intensity),...
%     'intensity is not a struct .');
% assert(ismember('width', fieldnames(intensity), ...
%     'intensity does not contain all of the required parameters');
% assert(isstruct(coord),...
%     'coord is not a struct .');
% assert(ismember('x', fieldnames(coord) && ...
%     ismember('', fieldnames(coord) && ...
%     ismember('', fieldnames(coord) && ...
%     ismember('', fieldnames(coord) && ...
%     ismember('', fieldnames(coord) && ...
%     ismember('', fieldnames(coord) && ...
%     ismember('', fieldnames(coord) && ...
% 
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
