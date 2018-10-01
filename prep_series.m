function [] = prep_series(result_file, image_path, image_names, ctrl_xw, ...
                  ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, fit_npts, ...
                  hue_lim, value_lim, entropy_lim, entropy_len, ...
                  morph_open_rad, morph_erode_rad, eql_len, xw, yw, mask_manual)
% function [] = prep_series(result_file, image_path, image_names, ctrl_xw, ...
%                   ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, fit_npts, ...
%                   hue_lim, value_lim, entropy_lim, entropy_len, ...
%                   morph_open_rad, morph_erode_rad, eql_len, xw, yw, mask_manual)
% 
% Create PIV input file for a given image series. Reads in the images, rectifies
% and crops, masks, corrects illuination, and saves the results and metadata in
% a netCDF file.
%
% Arguments:
% 
%   result_file = String, filename of the netCDF input file to be created. 
%
%   image_path = String, path to folder containing the images.
%
%   image_names = Cell array of strings, cells must contain filenames for
%       successive images in the experiment image series. 
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
%
% % Keith Ma

% check for sane arguments (pass-through arguments are checked in subroutines)
assert(nargin == 20);
validateattributes(result_file, {'char'}, {'vector'});
validateattributes(image_path, {'char'}, {'vector'});
validateattributes(image_names, {'cell'}, {'vector'});

% get some size parameters
nx = numel(xw);
ny = numel(yw);
num_image = numel(image_names);

% check that all images exist and have the expected size and type
info = imfinfo([image_path filesep image_names{1}]);
raw_nrow = info.Height;
raw_ncol = info.Width;
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
ncid = netcdf.create(result_file, 'NETCDF4');

% add global attributes 
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'yalebox git_hash', util_git_hash());
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop ctrl_xw', ctrl_xw);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop ctrl_yw', ctrl_yw);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop ctrl_xp', ctrl_xp);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop ctrl_yp', ctrl_yp);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop crop_xw', crop_xw);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop crop_yw', crop_yw);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_rectify_and_crop fit_npts', fit_npts);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_mask_auto hue_lim', hue_lim);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_mask_auto value_lim', value_lim);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_mask_auto entropy_lim', entropy_lim);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_mask_auto entropy_len', entropy_len);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_mask_auto morph_open_rad', morph_open_rad);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_mask_auto morph_erode_rad', morph_erode_rad);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'prep_intensity eql_len', eql_len);
     
% create dimensions
x_dimid = netcdf.defDim(ncid, 'x', nx);
y_dimid = netcdf.defDim(ncid, 'y', ny);
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

raw_varid = netcdf.defVar(ncid, 'img_rgb', 'NC_UBYTE', [y_dimid, x_dimid, rgb_dimid, step_dimid]);
netcdf.putAtt(ncid, raw_varid, 'long_name', 'rectified rgb image');
netcdf.putAtt(ncid, raw_varid, 'units', '24-bit color');
netcdf.defVarDeflate(ncid, raw_varid, true, true, 1);
netcdf.defVarChunking(ncid, raw_varid, 'CHUNKED', [ny, nx, 3, 1]);

img_varid = netcdf.defVar(ncid, 'img', 'NC_FLOAT', [y_dimid, x_dimid, step_dimid]);
netcdf.putAtt(ncid, img_varid, 'long_name', 'rectified normalized grayscale image');
netcdf.putAtt(ncid, img_varid, 'units', '1');
netcdf.defVarDeflate(ncid, img_varid, true, true, 1);
netcdf.defVarChunking(ncid, img_varid, 'CHUNKED', [ny, nx, 1]);

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
ncid = netcdf.open(result_file, 'WRITE');
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
     
    % update user (always verbose)
    fprintf('\n%s: %s\n', mfilename, this_file);
    
    % rectify and crop
    raw = prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, ...
              crop_yw, raw, fit_npts, false, true);
     
    % compute automatic mask
    mask_auto = prep_mask_auto(raw, hue_lim, value_lim, entropy_lim, ...
                    entropy_len, morph_open_rad, morph_erode_rad, false, true);
    
    % equalize intensity
    img = prep_intensity(raw, mask_manual & mask_auto, eql_len, false, true);
     
    % save results   
    ncid = netcdf.open(result_file, 'WRITE');
    netcdf.putVar(ncid, maska_varid, [0, 0, i-1], [ny, nx, 1], uint8(mask_auto));
    netcdf.putVar(ncid, img_varid, [0, 0, i-1], [ny, nx, 1], img);
    netcdf.putVar(ncid, raw_varid, [0, 0, 0, i-1], [ny, nx, 3, 1], raw);
    netcdf.close(ncid);
     
end
