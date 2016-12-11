function [] = post_series(piv_netcdf, post_netcdf, pro_bbox, ...
                          retro_bbox, pad_method, notes)
% 
% Run post-processing analyses on PIV data, and saves the results and
% metadata in a netCDF file.
%
% Arguments:
%   piv_netcdf: string, file name of netCDF containing PIV results as
%       produced by piv_series.m
%   post_netcdf: string, file name of netCDF to create for post-processing
%       results data and metadata. Overwriting is not allowed to prevent
%       unhappy accidents
%   pro_bbox, retro_bbox: 4-element position vectors [xmin, ymin, width,
%       height] for the bounding boxes used to estimate pro- and retro-
%       plate displacements from PIV data.
%   pad_method: String, ...
%   notes: String, notes to be included in output netCDF as a global
%       attribute. default = ''
%
% Output netCDF is self-documenting
%
% % Keith Ma

% set defaults
if nargin < 6
    notes = '';
end

% check (immediate) input arguments
validateattributes(piv_netcdf, {'char'}, {'nonempty'}, ...
    mfilename, 'piv_netcdf');
assert(exist(piv_netcdf, 'file') == 2, ...
    sprintf('input file %s does not exist', piv_netcdf));
validateattributes(post_netcdf, {'char'}, {'nonempty'}, ...
    mfilename, 'post_netcdf');
assert(exist(post_netcdf, 'file') ~= 2, ...
    sprintf('output file %s already exists', post_netcdf));
validateattributes(notes, {'char'}, {}, mfilename, 'notes');

% read in PIV data from netCDF
% % stored
xx   = double(ncread(piv_netcdf, 'x')); 
yy   = double(ncread(piv_netcdf, 'y')); 
step = double(ncread(piv_netcdf, 'step'));
uu   = double(ncread(piv_netcdf, 'u'));
vv   = double(ncread(piv_netcdf, 'v'));
roi  = logical(ncread(piv_netcdf, 'roi'));
% % derived
mm = sqrt(uu.^2+vv.^2);
num_steps = length(step);
num_x = length(xx);
num_y = length(yy);

% create output file
cmode = bitor(netcdf.getConstant('NETCDF4'), netcdf.getConstant('NOCLOBBER'));
ncid = netcdf.create(post_netcdf, cmode);

% add global attributes 
global_id = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid, global_id, 'yalebox commit hash', util_git_hash());
netcdf.putAtt(ncid, global_id, 'piv_netcdf MD5 hash', util_md5_hash(piv_netcdf));
netcdf.putAtt(ncid, global_id, 'pro_bbox', pro_bbox);
netcdf.putAtt(ncid, global_id, 'retro_bbox', retro_bbox);
netcdf.putAtt(ncid, global_id, 'pad_method', pad_method);
netcdf.putAtt(ncid, global_id, 'notes', notes);
% netcdf.putAtt(ncid, global_id, '', );

% define dimensions
x_dimid = netcdf.defDim(ncid, 'x', num_x);
y_dimid = netcdf.defDim(ncid, 'y', num_y);
s_dimid = netcdf.defDim(ncid, 'step', num_steps);

% define variables and thier attributes, compression, and chunking
dim_3d = [y_dimid, x_dimid, s_dimid];
chunk_3d = [num_y, num_x, 1];

x_varid = netcdf.defVar(ncid, 'x', 'NC_FLOAT', x_dimid);
netcdf.putAtt(ncid, x_varid, 'long_name', 'horizontal position');
netcdf.putAtt(ncid, x_varid, 'units', 'meters');

y_varid = netcdf.defVar(ncid, 'y', 'NC_FLOAT', y_dimid);
netcdf.putAtt(ncid, y_varid, 'long_name', 'vertical position');
netcdf.putAtt(ncid, y_varid, 'units', 'meters');

s_varid = netcdf.defVar(ncid, 'step', 'NC_FLOAT', s_dimid);
netcdf.putAtt(ncid, s_varid, 'long_name', 'step number');
netcdf.putAtt(ncid, s_varid, 'units', '1');

u_pro_varid = netcdf.defVar(ncid, 'u_pro', 'NC_FLOAT', s_dimid);
netcdf.putAtt(ncid, u_pro_varid, 'long_name', 'median proside section displacement vector, x-component');
netcdf.putAtt(ncid, u_pro_varid, 'units', 'meters/step');

v_pro_varid = netcdf.defVar(ncid, 'v_pro', 'NC_FLOAT', s_dimid);
netcdf.putAtt(ncid, v_pro_varid, 'long_name', 'median proside section displacement vector, y-component');
netcdf.putAtt(ncid, v_pro_varid, 'units', 'meters/step');

u_retro_varid = netcdf.defVar(ncid, 'u_retro', 'NC_FLOAT', s_dimid);
netcdf.putAtt(ncid, u_retro_varid, 'long_name', 'median retroside section displacement vector, x-component');
netcdf.putAtt(ncid, u_retro_varid, 'units', 'meters/step');

v_retro_varid = netcdf.defVar(ncid, 'v_retro', 'NC_FLOAT', s_dimid);
netcdf.putAtt(ncid, v_retro_varid, 'long_name', 'median retroside section displacement vector, y-component');
netcdf.putAtt(ncid, v_retro_varid, 'units', 'meters/step');

L11_varid = netcdf.defVar(ncid, 'L11', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, L11_varid, 'CHUNKED', chunk_3d);
netcdf.putAtt(ncid, L11_varid, 'long_name', 'du/dx, deformation gradient tensor element (1,1)');
netcdf.putAtt(ncid, L11_varid, 'units', '1');

L12_varid = netcdf.defVar(ncid, 'L12', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, L12_varid, 'CHUNKED', chunk_3d);
netcdf.putAtt(ncid, L12_varid, 'long_name', 'du/dy, deformation gradient tensor element (1,2)');
netcdf.putAtt(ncid, L12_varid, 'units', '1');

L21_varid = netcdf.defVar(ncid, 'L21', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, L21_varid, 'CHUNKED', chunk_3d);
netcdf.putAtt(ncid, L21_varid, 'long_name', 'dv/dx, deformation gradient tensor element (2,1)');
netcdf.putAtt(ncid, L21_varid, 'units', '1');

L22_varid = netcdf.defVar(ncid, 'L22', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, L22_varid, 'CHUNKED', chunk_3d);
netcdf.putAtt(ncid, L22_varid, 'long_name', 'dv/dy, deformation gradient tensor element (2,2)');
netcdf.putAtt(ncid, L22_varid, 'units', '1');

% finish netcdf creation
netcdf.endDef(ncid);
netcdf.close(ncid);    

% populate dimension values
ncid = netcdf.open(post_netcdf, 'WRITE');
netcdf.putVar(ncid, x_varid, xx);
netcdf.putVar(ncid, y_varid, yy);
netcdf.putVar(ncid, s_varid, step);
netcdf.close(ncid);

% run analyses for each timestep
for ii = 1:num_steps
    us = uu(:,:,ii);
    vs = vv(:,:,ii);
    rois = roi(:,:,ii);
    uv_pro   = post_displ_rect(xx, yy, us, vs, pro_bbox);
    uv_retro = post_displ_rect(xx, yy, us, vs, retro_bbox);
    [L, F] = post_strain(xx, yy, us, vs, rois, pad_method, true);
    
    % save results
    ncid = netcdf.open(post_netcdf, 'WRITE');
    netcdf.putVar(ncid, u_pro_varid, ii-1, 1, uv_pro(1));
    netcdf.putVar(ncid, v_pro_varid, ii-1, 1, uv_pro(2));    
    netcdf.putVar(ncid, u_retro_varid, ii-1, 1, uv_retro(1));
    netcdf.putVar(ncid, v_retro_varid, ii-1, 1, uv_retro(2));
    netcdf.putVar(ncid, L11_varid, [0, 0, ii-1], [num_y, num_x, 1], L(:,:,1));
    netcdf.putVar(ncid, L21_varid, [0, 0, ii-1], [num_y, num_x, 1], L(:,:,2));
    netcdf.putVar(ncid, L12_varid, [0, 0, ii-1], [num_y, num_x, 1], L(:,:,3));
    netcdf.putVar(ncid, L22_varid, [0, 0, ii-1], [num_y, num_x, 1], L(:,:,4));
    netcdf.close(ncid);
    break
end

% finalize

% debug
keyboard