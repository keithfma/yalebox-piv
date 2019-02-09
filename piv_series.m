function [] = piv_series(...
    output_file, input_file, step_range, gap, samp_len, samp_spc, ...
    intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
    min_frac_data, min_frac_overlap, verbose, version)
%
% function [] = piv_series(...
%     output_file, input_file, step_range, gap, samp_len, samp_spc, ...
%     intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
%     min_frac_data, min_frac_overlap, verbose, version)
%
% Run PIV analysis for a given input series. Input image data  expected to
% be a MAT file as created by prep_series(). Results are saved in a new
% MAT file which includes all relevant metadata.
%
% Arguments:
%   output_file: String, name of the output MAT file containing PIV results
%   input_file: String, name of the input MAT file containing pre-processed
%       image data
%   step_range: 2-element vector, [initial, final] steps of input images series
%       to include in PIV analysis. Either element can be set to NaN to use the
%       full range.
%   gap : Scalar, integer, gap between analyzed images in steps, for
%       example, for an initial image at step 3, setting gap -> 1 would use
%       a final image at step 4 and yield PIV results at step 3.5, or
%       alternatively, setting gap -> 2 would use a final image at step 5
%       and yield PIV results at step 4.
%   samp_len: see piv() help for details
%   samp_spc: see piv() help for details
%   intr_len: see piv() help for details
%   num_pass: see piv() help for details
%   valid_radius: see piv() help for details
%   valid_max: see piv() help for details
%   valid_eps: see piv() help for details
%   min_frac_data: see piv() help for details
%   min_frac_overlap: see piv() help for details
%   verbose: Scalar, logical, display verbose messages for this function and its
%       children (1) or don't (0)
%   version: String, optional, version hash for yalebox-piv, required as a
%       work-around for path problems on the SCC, attempt to read
%       automatically if not provided.
% %

update_path('prep', 'util');

% check for sane arguments (pass-through arguments are checked in subroutines)
validateattributes(output_file, {'char'}, {'vector'});
validateattributes(input_file, {'char'}, {'vector'});
validateattributes(step_range, {'numeric'}, {'vector', 'numel', 2});
validateattributes(gap, {'numeric'}, {'scalar', 'positive', 'integer'});
[~, ~, result_file_ext] = fileparts(result_file);
assert(strcmp('.mat', result_file_ext), 'Output file must be .mat');

% read input dimension values
x_img = double(ncread(input_file,'x'));
y_img = double(ncread(input_file,'y'));
step_img = double(ncread(input_file, 'step'));

% compute indices for image pairs
if isnan(step_range(1))
    min_idx = 1;
else
    min_idx = find(step_img >= step_range(1), 1, 'first');
end
if isnan(step_range(2))
    max_idx = length(step_img);
else
    max_idx = find(step_img <= step_range(2), 1, 'last');
end
ini_idx = 1:length(step_img);
fin_idx = ini_idx + gap;
in_range = ini_idx >= min_idx & ini_idx <= max_idx ...
         & fin_idx >= min_idx & fin_idx <= max_idx;
ini_idx = ini_idx(in_range);
fin_idx = fin_idx(in_range);

% compute size of gridded output dimensions
[r_grd, c_grd] = piv_sample_grid(samp_spc, length(y_img), length(x_img));
step = 0.5*(step_img(ini_idx) + step_img(fin_idx));
num_x_grd = size(c_grd, 2);
num_y_grd = size(r_grd, 1);
num_step = length(step);
max_num_pts = num_x_grd*num_y_grd;

% create netcdf file
ncid = netcdf.create(output_file, 'NETCDF4');

% add global attributes
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv_series step_range', step_range);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv_series gap', gap);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv samp_len', samp_len);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv samp_spc', samp_spc);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv intr_len', intr_len);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv num_pass', num_pass);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv valid_radius', valid_radius);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv valid_max', valid_max);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv valid_eps', valid_eps);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv min_frac_data', min_frac_data);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv min_frac_overlap', min_frac_overlap);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'yalebox-piv version', version);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'image file MD5 hash', util_md5_hash(input_file));

% create dimensions
x_grd_dimid = netcdf.defDim(ncid, 'x_grd', num_x_grd);
y_grd_dimid = netcdf.defDim(ncid, 'y_grd', num_y_grd);
s_dimid = netcdf.defDim(ncid, 'step', num_step);
p_dimid = netcdf.defDim(ncid, 'pts', netcdf.getConstant('NC_UNLIMITED'));

% define gridded dimensions and variables
grd_dims = [y_grd_dimid, x_grd_dimid, s_dimid];
grd_chunks = [num_y_grd, num_x_grd, 1];

x_grd_varid = netcdf.defVar(ncid, 'x_grd', 'NC_DOUBLE', x_grd_dimid);
netcdf.putAtt(ncid, x_grd_varid, 'long_name', 'horizontal position, regular grid');
netcdf.putAtt(ncid, x_grd_varid, 'units', 'meters');

y_grd_varid = netcdf.defVar(ncid, 'y_grd', 'NC_DOUBLE', y_grd_dimid);
netcdf.putAtt(ncid, y_grd_varid, 'long_name', 'vertical position, regular grid');
netcdf.putAtt(ncid, y_grd_varid, 'units', 'meters');

s_varid = netcdf.defVar(ncid, 'step', 'NC_DOUBLE', s_dimid);
netcdf.putAtt(ncid, s_varid, 'long_name', 'step number');
netcdf.putAtt(ncid, s_varid, 'units', '1');

u_grd_varid = netcdf.defVar(ncid, 'u_grd', 'NC_DOUBLE', grd_dims);
netcdf.putAtt(ncid, u_grd_varid, 'long_name', ...
    'displacement vector, x-component, interpolated to regular grid');
netcdf.putAtt(ncid, u_grd_varid, 'units', 'meters/step');
netcdf.defVarDeflate(ncid, u_grd_varid, true, true, 1);
netcdf.defVarChunking(ncid, u_grd_varid, 'CHUNKED', grd_chunks);

v_grd_varid = netcdf.defVar(ncid, 'v_grd', 'NC_DOUBLE', grd_dims);
netcdf.putAtt(ncid, v_grd_varid, 'long_name', ...
    'displacement vector, y-component, interpolated to regular grid');
netcdf.putAtt(ncid, v_grd_varid, 'units', 'meters/step');
netcdf.defVarDeflate(ncid, v_grd_varid, true, true, 1);
netcdf.defVarChunking(ncid, v_grd_varid, 'CHUNKED', grd_chunks);

r_grd_varid = netcdf.defVar(ncid, 'roi_grd', 'NC_BYTE', grd_dims);
netcdf.putAtt(ncid, r_grd_varid, 'long_name', ...
    'displacement vector mask, regular grid');
netcdf.putAtt(ncid, r_grd_varid, 'units', 'boolean');
netcdf.defVarDeflate(ncid, r_grd_varid, true, true, 1);
netcdf.defVarChunking(ncid, r_grd_varid, 'CHUNKED', grd_chunks);

% define scattered-point variables, thier attributes, compression, and chunking
pts_dims = [p_dimid, s_dimid];
pts_chunks = [max_num_pts, 1];
pts_fill = NaN;

x_pts_varid = netcdf.defVar(ncid, 'x_pts', 'NC_DOUBLE', pts_dims);
netcdf.putAtt(ncid, x_pts_varid, 'long_name', ...
    'horizontal position for raw measurements at scattered points');
netcdf.putAtt(ncid, x_pts_varid, 'units', 'meters');
netcdf.defVarDeflate(ncid, x_pts_varid, true, true, 1);
netcdf.defVarChunking(ncid, x_pts_varid, 'CHUNKED', pts_chunks);
netcdf.defVarFill(ncid, x_pts_varid, false, pts_fill); 

y_pts_varid = netcdf.defVar(ncid, 'y_pts', 'NC_DOUBLE', pts_dims);
netcdf.putAtt(ncid, y_pts_varid, 'long_name', ...
    'horizontal position for raw measurements at scattered points');
netcdf.putAtt(ncid, y_pts_varid, 'units', 'meters');
netcdf.defVarDeflate(ncid, y_pts_varid, true, true, 1);
netcdf.defVarChunking(ncid, y_pts_varid, 'CHUNKED', pts_chunks);
netcdf.defVarFill(ncid, y_pts_varid, false, pts_fill);

u_pts_varid = netcdf.defVar(ncid, 'u_pts', 'NC_DOUBLE', pts_dims);
netcdf.putAtt(ncid, u_pts_varid, 'long_name', ...
    'displacement vector, x-component, raw measurements at scattered points');
netcdf.putAtt(ncid, u_pts_varid, 'units', 'meters');
netcdf.defVarDeflate(ncid, u_pts_varid, true, true, 1);
netcdf.defVarChunking(ncid, u_pts_varid, 'CHUNKED', pts_chunks);
netcdf.defVarFill(ncid, u_pts_varid, false, pts_fill);

v_pts_varid = netcdf.defVar(ncid, 'v_pts', 'NC_DOUBLE', pts_dims);
netcdf.putAtt(ncid, v_pts_varid, 'long_name', ...
    'displacement vector, y-component, raw measurements at scattered points');
netcdf.putAtt(ncid, v_pts_varid, 'units', 'meters');
netcdf.defVarDeflate(ncid, v_pts_varid, true, true, 1);
netcdf.defVarChunking(ncid, v_pts_varid, 'CHUNKED', pts_chunks);
netcdf.defVarFill(ncid, v_pts_varid, false, pts_fill);

% finish netcdf creation
netcdf.endDef(ncid);
netcdf.close(ncid);

% read constants
roi_const = logical(ncread(input_file, 'mask_manual', [1, 1], [inf, inf]));

% analyse all steps
for ii = 1:num_step
    
    if verbose
        fprintf('\n%s: begin step = %.1f\n', mfilename, step(ii));
        fprintf('%s: ini_step = %.1f\n', mfilename, step_img(ini_idx(ii)));
        fprintf('%s: fin_step = %.1f\n', mfilename, step_img(fin_idx(ii)));
    end
    
    % update image and roi pair
    img0 = double(ncread(input_file, 'img', [1, 1, ini_idx(ii)], [inf, inf, 1]));
    roi0 = logical(ncread(input_file, 'mask_auto', [1, 1, ini_idx(ii)], [inf, inf, 1])) & roi_const;
    
    img1 = double(ncread(input_file, 'img', [1, 1, fin_idx(ii)], [inf, inf, 1]));
    roi1 = logical(ncread(input_file, 'mask_auto', [1, 1, fin_idx(ii)], [inf, inf, 1])) & roi_const;
    
    % perform piv analysis
    result = piv(img0, img1, roi0, roi1, x_img, y_img, samp_len, samp_spc, ...
        intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
        min_frac_data, min_frac_overlap, verbose); 
    
    % write results to output file
    ncid = netcdf.open(output_file, 'WRITE');
    
    grd_start = [0, 0, ii-1];
    grd_count = [num_y_grd, num_x_grd, 1];
    
    netcdf.putVar(ncid, u_grd_varid, grd_start, grd_count, result.u_grd);
    netcdf.putVar(ncid, v_grd_varid, grd_start, grd_count, result.v_grd); 
    netcdf.putVar(ncid, r_grd_varid, grd_start, grd_count, int8(result.roi_grd));
    
    pts_start = [0, ii-1];
    pts_count = [length(result.x_pts), 1];
    
    netcdf.putVar(ncid, x_pts_varid, pts_start, pts_count, result.x_pts);
    netcdf.putVar(ncid, y_pts_varid, pts_start, pts_count, result.y_pts);
    netcdf.putVar(ncid, u_pts_varid, pts_start, pts_count, result.u_pts);
    netcdf.putVar(ncid, v_pts_varid, pts_start, pts_count, result.v_pts);
    netcdf.close(ncid);
    
    if verbose
        fprintf('%s: end step = %.1f\n', mfilename,  step(ii));
    end
    
end

% populate dimension values
ncid = netcdf.open(output_file, 'WRITE');
netcdf.putVar(ncid, x_grd_varid, result.x_grd(1,:));
netcdf.putVar(ncid, y_grd_varid, result.y_grd(:,1));
netcdf.putVar(ncid, s_varid, step);
netcdf.close(ncid);