% function [] = prep_series(output_file, image_path, image_names, verbose)
% function [] = prep_series(output_file, image_path, image_names, x, y, scale, ...
%                   offset, mask_manual, hue_lim, val_lim, entr_lim, entr_win, ...
%                   morph_open_rad, morph_erode_rad, nwin, verbose)
% 
% Create PIV input file for a given image series. Reads in the images,
% performs masking and color correction, and saves the results and metadata
% in a netCDF file.
%
% Arguments:
% 
%   output_file = String, filename of the netCDF input file to be created. 
%
%   image_path = String, path to folder containing the images.
%
%   image_names = Cell array of strings, cells must contain filenames for
%       successive images in the experiment image series. 
%
%   x, y, scale, offset = Output arguments from prep_world_coord().
%
%   mask_manual = Output argument from prep_mask_manual()
%
%   hue_lim, val_lim, entr_lim, entr_win, morph_open_rad, morph_erode_rad = Input 
%       arguments for prep_mask_auto()
%
%   nwin = Input argument from prep_intensity()
%
%   verbose = Logical flag, display verbose messages (1) or don't (0)
%
%   PIV input netCDF format:
%       dimensions: x, y, step
%       variables: mask_auto, mask_manual, x, y, step, intensity
%       attributes: all preprocessing parameters
%
% %

% DEBUG
output_file = 'junk.nc';
image_names = {'K24_050.jpg', 'K24_100.jpg', 'K24_150.jpg'};
image_path = '/home/kfm/Documents/dissertation/yalebox-exp-erosion/data/k24/image/clean';
verbose = true;

% % check for sane arguments (pass-through arguments are checked in subroutines)
% narginchk(16, 16); 
% validateattributes(output_file, {'char'}, {'vector'});
% validateattributes(image_path, {'char'}, {'vector'});
% validateattributes(image_names, {'cell'}, {'vector'});
% 

% get some size parameters
nx = numel(xw);
ny = numel(yw);
info = imfinfo([image_path filesep image_names{1}]);
raw_nrow = info.Height;
raw_ncol = info.Width;
num_image = numel(image_names);

% check that all images exist and have the expected size and type
for i = 1:num_image
    try
        this_file = [image_path filesep image_names{i}];
        info = imfinfo(this_file);
        assert(info.Width == raw_ncol && info.Height == raw_nrow, ...
            sprintf('incorrect dimensions in image %s', this_file));
        assert(info.BitDepth == 24, ...
            sprintf('incorrect bit depth in image %s', this_file));
    catch
        error('unable to read image %i: %s', i, this_file);
    end
end

% create netcdf file
ncid = netcdf.create(output_file, 'NETCDF4');

% add global attributes 
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop ctrl_xw', ctrl_xw);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop ctrl_yw', ctrl_yw);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop ctrl_xp', ctrl_xp);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop ctrl_yp', ctrl_yp);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop crop_xw', crop_xw);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop crop_yw', crop_yw);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_mask_auto entropy_len', entropy_len);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_mask_auto num_cluster', num_cluster);
for ii = 1:num_cluster
    netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), sprintf('prep_mask_auto cluster_center_%d', ii), cluster_center(ii,:));
end
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_intensity eql_len', eql_len);
     
% create dimensions
x_dimid = netcdf.defDim(ncid, 'x', nx);
y_dimid = netcdf.defDim(ncid, 'y', ny);
row_dimid = netcdf.defDim(ncid, 'row', raw_nrow);
col_dimid = netcdf.defDim(ncid, 'col', raw_ncol);
step_dimid = netcdf.defDim(ncid, 'step', numel(image_names));
rgb_dimid = netcdf.defDim(ncid, 'rgb', 3);
 
% define variables 
x_varid = netcdf.defVar(ncid, 'x', 'NC_FLOAT', x_dimid);
netcdf.putAtt(ncid, x_varid, 'long_name', 'horizontal position');
netcdf.putAtt(ncid, x_varid, 'units', 'meters');
 
y_varid = netcdf.defVar(ncid, 'y', 'NC_FLOAT', y_dimid);
netcdf.putAtt(ncid, y_varid, 'long_name', 'vertical position');
netcdf.putAtt(ncid, y_varid, 'units', 'meters');
 
step_varid = netcdf.defVar(ncid, 'step', 'NC_SHORT', step_dimid);
netcdf.putAtt(ncid, step_varid, 'long_name', 'step number');
netcdf.putAtt(ncid, step_varid, 'units', '1');
 
rgb_varid = netcdf.defVar(ncid, 'rgb', 'NC_CHAR', rgb_dimid);
netcdf.putAtt(ncid, rgb_varid, 'long_name', 'color band');
netcdf.putAtt(ncid, rgb_varid, 'units', 'char');

raw_varid = netcdf.defVar(ncid, 'raw', 'NC_UBYTE', [row_dimid, col_dimid, rgb_dimid, step_dimid]);
netcdf.putAtt(ncid, raw_varid, 'long_name', 'original image');
netcdf.putAtt(ncid, raw_varid, 'units', '24-bit color');
netcdf.defVarDeflate(ncid, raw_varid, true, true, 1);
netcdf.defVarChunking(ncid, raw_varid, 'CHUNKED', [raw_nrow, raw_ncol, 3, 1]);

image_varid = netcdf.defVar(ncid, 'image', 'NC_FLOAT', [y_dimid, x_dimid, step_dimid]);
netcdf.putAtt(ncid, image_varid, 'long_name', 'normalized sand brightness');
netcdf.putAtt(ncid, image_varid, 'units', '1');
netcdf.defVarDeflate(ncid, image_varid, true, true, 1);
netcdf.defVarChunking(ncid, image_varid, 'CHUNKED', [ny, nx, 1]);

maska_varid = netcdf.defVar(ncid, 'mask_auto', 'NC_BYTE', [y_dimid, x_dimid, step_dimid]);
netcdf.putAtt(ncid, maska_varid, 'long_name', 'sand mask, automatic');
netcdf.putAtt(ncid, maska_varid, 'units', 'boolean');
netcdf.defVarDeflate(ncid, maska_varid, true, true, 1);
netcdf.defVarChunking(ncid, maska_varid, 'CHUNKED', [ny, nx, 1]);
 
maskm_varid = netcdf.defVar(ncid, 'mask_manual', 'NC_BYTE', [y_dimid, x_dimid]);
netcdf.putAtt(ncid, maskm_varid, 'long_name', 'sand mask, manual');
netcdf.putAtt(ncid, maskm_varid, 'units', 'boolean');
netcdf.defVarDeflate(ncid, maskm_varid, true, true, 1);

% finish netcdf creation
netcdf.endDef(ncid);
netcdf.close(ncid);

% populate constant variables
ncid = netcdf.open(output_file, 'WRITE');
netcdf.putVar(ncid, x_varid, xw);
netcdf.putVar(ncid, y_varid, yw);
netcdf.putVar(ncid, rgb_varid, 'rgb');
netcdf.putVar(ncid, step_varid, 0:num_image-1);
netcdf.putVar(ncid, maskm_varid, uint8(mask_manual));
netcdf.close(ncid);
 
% loop over all images
for i = 1:num_image
    
    % read in original image
    this_file = [image_path filesep image_names{i}];
    raw = imread(this_file);
     
    if verbose
        fprintf('%s: %s\n', mfilename, this_file);
    end
    
    % rectify and crop
    rgb = prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, ...
              crop_yw, raw, false, true);
     
    % compute automatic mask
    mask_auto = prep_mask_auto(rgb, entropy_len, cluster_center, false, true);

    % equalize intensity
    hsv = rgb2hsv(rgb);
    value = hsv(:,:,3);
    eql = prep_intensity(value, mask_manual & mask_auto, eql_len, false, true);
     
    % save results   
    ncid = netcdf.open(output_file, 'WRITE');
    netcdf.putVar(ncid, maska_varid, [0, 0, i-1], [ny, nx, 1], ...
        uint8(mask_auto));
    netcdf.putVar(ncid, image_varid, [0, 0, i-1], [ny, nx, 1], ...
        eql);
    netcdf.putVar(ncid, raw_varid, [0, 0, 0, i-1], [raw_nrow, raw_ncol, 3, 1], ...
        raw);
    netcdf.close(ncid);
     
end
