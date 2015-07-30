function [] = yalebox_prep_create_input(input_file, image_files, param)
% 
% Create PIV input file for a given image series. Reads in the images,
% performs masking and color correction, and saves the results and metadata
% in a netCDF file.
%
% All input parameters for the subroutines are included in the param
% struct, and must be explicitly specified. Typically, one would first
% experiment on a subset of images from an experiment to determine
% successful settings, then construct the param mask and run this function
% to process the entire series.
%
% Arguments:
% 
% input_file = 
% images_files = 
% param = 
%
% PIV input file format:
%
%
% Keith Ma, July 2015

% check input argument types
assert(isa(input_file, 'char'), ...
    'input_file is not a string.');
assert(isa(image_files, 'cell'), ...
    'image_files is not a cell array');
assert(isstruct(param),...
    'param is not a struct.');

% check that all required parameters exist in the param struct

% check that all of the image files exist

% check that image and mask dimensions match the coordinate vectors


% for i = 1:numel(image_files)
%     assert(exist(image_files{i}, 'file') == 2, ...
%         sprintf('image_files{%i} does not exist', i));
%    
%     info = iminfo(image_files{i});
%     if i>1
%         assert(info.Width == info_prev.Width && info.Height == info_prev.Height, ...
%             sprintf('image_files{%i} is not the same size as the previous image', i));
%     end
%     info_prev = info;
% end
%         
        
% CHECK THE CONTENTS OF THE PARAM MATRIX

% init netcdf file
% ADD PARAMETERS
% ADD DIMENSIONS AND VARIABLES

% loop over all images

% read in original image
