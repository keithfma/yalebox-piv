function [] = yalebox_prep_create_input(input_file, image_path, image_names, x, y, x_scale, y_scale, x_offset, y_offset, mask_manual, hue_lim, value_lim, entropy_lim, median_window, entropy_window, morph_radius, histeq_width)
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
%   morph_radius = Select input arguments from yalebox_prep_mask_auto()
% 
% histeq_width = Select input argument from yalebox_prep_intensity()
%
% PIV input netCDF format:
%   groups: preprocess
%   dimensions: x, y, step
%   variables: preprocess/mask_auto, preprocess/mask_manual, x, y, step, intensity
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

% check that coordinate vectors match image dimensions
image_w = size(mask_manual, 2);
image_h = size(mask_manual, 1);
assert(numel(x) == image_w && numel(y) == image_h, ...
    'coordinate vectors and image are not the same size');

% check that all images exist and have the expected size and type
nimage = numel(image_names);
for i = 1:nimage
    try
        this_file = [image_path filesep image_names{i}];
        info = imfinfo(this_file);
        assert(info.Width == image_w && info.Height == image_h, ...
            sprintf('incorrect dimensions in image %s', this_file));
        assert(info.BitDepth == 24, ...
            sprintf('incorrect bit depth in image %s', this_file));
    catch 
        error('Unable to read image %i: %s', i, this_file);
    end
end

% create netcdf file
ncid = netcdf.create(input_file, 'NETCDF4');

% create dimensions
x_dimid = netcdf.defDim(ncid, 'x', numel(x));
y_dimid = netcdf.defDim(ncid, 'y', numel(y));
s_dimid = netcdf.defDim(ncid, 'step', numel(image_names));

% create groups with attributes
p_grpid = netcdf.defGrp(ncid, 'preprocess');
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_world_coord x_scale', x_scale);
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_world_coord y_scale', y_scale);
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_world_coord x_offset', x_offset);
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_world_coord y_offset', y_offset);
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_mask_auto hue_lim', hue_lim);
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_mask_auto value_lim', value_lim);
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_mask_auto entropy_lim', entropy_lim);
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_mask_auto median_window', median_window);
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_mask_auto entropy_window', entropy_window);
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_mask_auto morph_radius', morph_radius);
netcdf.putAtt(p_grpid, netcdf.getConstant('GLOBAL'),...
    'yalebox_prep_intensity histeq_width', histeq_width);

% mask_manual

% define variables with attributes
D3 = [x_dimid, y_dimid, s_dimid];
D2 = [x_dimid, y_dimid];

x_varid = netcdf.defVar(ncid, 'x', 'NC_DOUBLE', x_dimid);
netcdf.putAtt(ncid, x_varid, 'long_name', 'horizontal position');
netcdf.putAtt(ncid, x_varid, 'units', 'meters');

y_varid = netcdf.defVar(ncid, 'y', 'NC_DOUBLE', y_dimid);
netcdf.putAtt(ncid, y_varid, 'long_name', 'vertical position');
netcdf.putAtt(ncid, y_varid, 'units', 'meters');

s_varid = netcdf.defVar(ncid, 'step', 'NC_DOUBLE', s_dimid);
netcdf.putAtt(ncid, s_varid, 'long_name', 'step number');
netcdf.putAtt(ncid, s_varid, 'units', '1');

i_varid = netcdf.defVar(ncid, 'intensity', 'NC_DOUBLE', D3);
netcdf.putAtt(ncid, i_varid, 'long_name', 'normalized sand brightness');
netcdf.putAtt(ncid, i_varid, 'units', '1');

ma_varid = netcdf.defVar(p_grpid, 'mask_auto', 'NC_BYTE', D3);
netcdf.putAtt(p_grpid, ma_varid, 'long_name', 'sand mask, automatic');
netcdf.putAtt(p_grpid, ma_varid, 'units', 'boolean');

mm_varid = netcdf.defVar(p_grpid, 'mask_manual', 'NC_BYTE', D2);
netcdf.putAtt(p_grpid, mm_varid, 'long_name', 'sand mask, manual');
netcdf.putAtt(p_grpid, mm_varid, 'units', 'boolean');

% finish netcdf creation
netcdf.endDef(ncid);
netcdf.close(ncid);

% populate constant variables
ncid = netcdf.open(input_file, 'WRITE');
netcdf.putVar(ncid, x_varid, x);
netcdf.putVar(ncid, y_varid, y);
netcdf.putVar(ncid, s_varid, 0:nimage-1);
netcdf.putVar(p_grpid, mm_varid, uint8(mask_manual'));
netcdf.close(ncid);

% loop over all images
for i = 1:nimage
    
    
    % read in original image
    this_file = [image_path filesep image_names{i}];
    rgb = imread(this_file);
    
    % compute automatic mask
    
    % merge (constant) manual mask
    
    % equalize intensity
    
    % debug
    mask_auto = i*ones(size(mask_manual));
    intensity = i*ones(size(mask_manual));
    
    % save results
    ncid = netcdf.open(input_file, 'WRITE');
    netcdf.putVar(p_grpid, ma_varid, [0, 0, i-1], [image_w, image_h, 1], uint8(mask_auto'));
    netcdf.putVar(ncid, i_varid, [0, 0, i-1], [image_w, image_h, 1], intensity');
    netcdf.close(ncid);
    
end
% end loop