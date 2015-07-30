function [] = yalebox_prep_create_input(input_file, image_files, param)
% 
% Create PIV input file for a given image series. Reads in the images,
% performs masking and color correction, and saves the results and metadata
% in a netCDF file.
%
% Arguments:
% 
% input_file = 
% images_files = 
% param = 
%
% Input file format:
%

% check for sane inputs
assert(isa(image_files, 'cell'), ...
    'image_files is not a cell array');

for i = 1:numel(image_files)
    assert(exist(image_files{i}, 'file') == 2, ...
        sprintf('image_files{%i} does not exist', i));
   
    info = iminfo(image_files{i});
    if i>1
        assert(info.Width == info_prev.Width && info.Height == info_prev.Height, ...
            sprintf('image_files{%i} is not the same size as the previous image', i));
    end
    info_prev = info;
end
        
        
% CHECK THE CONTENTS OF THE PARAM MATRIX

% init netcdf file
% ADD PARAMETERS
% ADD DIMENSIONS AND VARIABLES

% loop over all images

% read in original image
