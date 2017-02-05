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
cmp_lvl = 2; 

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

F11_varid = netcdf.defVar(ncid, 'F11', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, F11_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, F11_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, F11_varid, 'long_name', 'dx/dX, deformation gradient tensor element (1,1)');
netcdf.putAtt(ncid, F11_varid, 'units', '1');

F12_varid = netcdf.defVar(ncid, 'F12', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, F12_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, F12_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, F12_varid, 'long_name', 'dx/dY, deformation gradient tensor element (1,2)');
netcdf.putAtt(ncid, F12_varid, 'units', '1');

F21_varid = netcdf.defVar(ncid, 'F21', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, F21_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, F21_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, F21_varid, 'long_name', 'dy/dX, deformation gradient tensor element (2,1)');
netcdf.putAtt(ncid, F21_varid, 'units', '1');

F22_varid = netcdf.defVar(ncid, 'F22', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, F22_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, F22_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, F22_varid, 'long_name', 'dy/dY, deformation gradient tensor element (2,2)');
netcdf.putAtt(ncid, F22_varid, 'units', '1');

S1_varid = netcdf.defVar(ncid, 'S1', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, S1_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, S1_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, S1_varid, 'long_name', 'maximum principle stretch');
netcdf.putAtt(ncid, S1_varid, 'units', '1');

S1x_varid = netcdf.defVar(ncid, 'S1x', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, S1x_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, S1x_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, S1x_varid, 'long_name', 'maximum principle stretch unit direction vector, x-component');
netcdf.putAtt(ncid, S1x_varid, 'units', '1');

S1y_varid = netcdf.defVar(ncid, 'S1y', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, S1y_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, S1y_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, S1y_varid, 'long_name', 'maximum principle stretch unit direction vector, y-component');
netcdf.putAtt(ncid, S1y_varid, 'units', '1');

S2_varid = netcdf.defVar(ncid, 'S2', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, S2_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, S2_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, S2_varid, 'long_name', 'minimum principle stretch');
netcdf.putAtt(ncid, S2_varid, 'units', '1');

S2x_varid = netcdf.defVar(ncid, 'S2x', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, S2x_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, S2x_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, S2x_varid, 'long_name', 'minimum principle stretch unit direction vector, x-component');
netcdf.putAtt(ncid, S2x_varid, 'units', '1');

S2y_varid = netcdf.defVar(ncid, 'S2y', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, S2y_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, S2y_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, S2y_varid, 'long_name', 'minimum principle stretch unit direction vector, y-component');
netcdf.putAtt(ncid, S2y_varid, 'units', '1');

spin_varid = netcdf.defVar(ncid, 'spin', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, spin_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, spin_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, spin_varid, 'long_name', 'Local rotation angle, clockwise');
netcdf.putAtt(ncid, spin_varid, 'units', 'radians');

D1_varid = netcdf.defVar(ncid, 'D1', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, D1_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, D1_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, D1_varid, 'long_name', 'equivalent maximum principle stretch rate');
netcdf.putAtt(ncid, D1_varid, 'units', '1/step');

D2_varid = netcdf.defVar(ncid, 'D2', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, D2_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, D2_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, D2_varid, 'long_name', 'equivalent minimum principle stretch rate');
netcdf.putAtt(ncid, D2_varid, 'units', '1/step');

Dd_varid = netcdf.defVar(ncid, 'Dd', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, Dd_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, Dd_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, Dd_varid, 'long_name', 'deviatoric strain rate');
netcdf.putAtt(ncid, Dd_varid, 'units', '1/step');
netcdf.putAtt(ncid, Dd_varid, 'references', 'Brandon (1995), Ring and Brandon (1999)'); 

Dv_varid = netcdf.defVar(ncid, 'Dv', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, Dv_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, Dv_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, Dv_varid, 'long_name', 'volume strain rate');
netcdf.putAtt(ncid, Dv_varid, 'units', '1/step');
netcdf.putAtt(ncid, Dv_varid, 'references', 'Brandon (1995), Ring and Brandon (1999)'); 

Wk_star_varid = netcdf.defVar(ncid, 'Wk_star', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, Wk_star_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, Wk_star_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, Wk_star_varid, 'long_name', 'kinematic dilatancy number, relative to Dd');
netcdf.putAtt(ncid, Wk_star_varid, 'units', '1');
netcdf.putAtt(ncid, Wk_star_varid, 'references', 'Brandon (1995), Ring and Brandon (1999)'); 

Ak_star_varid = netcdf.defVar(ncid, 'Ak_star', 'NC_FLOAT', dim_3d);
netcdf.defVarChunking(ncid, Ak_star_varid, 'CHUNKED', chunk_3d);
netcdf.defVarDeflate(ncid, Ak_star_varid, true, true, cmp_lvl);
netcdf.putAtt(ncid, Ak_star_varid, 'long_name', 'kinematic vorticity number, relative to Dd');
netcdf.putAtt(ncid, Ak_star_varid, 'units', '1');
netcdf.putAtt(ncid, Ak_star_varid, 'references', 'Brandon (1995), Ring and Brandon (1999)'); 

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
    
    % get results
    us = uu(:,:,ii);
    vs = vv(:,:,ii);
    rois = roi(:,:,ii);
    uv_pro   = post_displ_rect(xx, yy, us, vs, pro_bbox);
    uv_retro = post_displ_rect(xx, yy, us, vs, retro_bbox);
    strain = post_strain(xx, yy, us, vs, rois, pad_method);
    
    % save results
    ncid = netcdf.open(post_netcdf, 'WRITE');
    netcdf.putVar(ncid, u_pro_varid, ii-1, 1, uv_pro(1));
    netcdf.putVar(ncid, v_pro_varid, ii-1, 1, uv_pro(2));    
    netcdf.putVar(ncid, u_retro_varid, ii-1, 1, uv_retro(1));
    netcdf.putVar(ncid, v_retro_varid, ii-1, 1, uv_retro(2));
    netcdf.putVar(ncid, F11_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.F11);
    netcdf.putVar(ncid, F21_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.F21);
    netcdf.putVar(ncid, F12_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.F12);
    netcdf.putVar(ncid, F22_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.F22);
    netcdf.putVar(ncid, S1_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.S1);
    netcdf.putVar(ncid, S1x_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.S1x);
    netcdf.putVar(ncid, S1y_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.S1y);
    netcdf.putVar(ncid, S2_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.S2);
    netcdf.putVar(ncid, S2x_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.S2x);
    netcdf.putVar(ncid, S2y_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.S2y);
    netcdf.putVar(ncid, spin_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.spin);
    netcdf.putVar(ncid, D1_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.D1);
    netcdf.putVar(ncid, D2_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.D2);
    netcdf.putVar(ncid, Dd_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.Dd);
    netcdf.putVar(ncid, Dv_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.Dv);
    netcdf.putVar(ncid, Wk_star_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.Wk_star);
    netcdf.putVar(ncid, Ak_star_varid, [0, 0, ii-1], [num_y, num_x, 1], strain.Ak_star);
    netcdf.close(ncid);
    
    % report progress
    fprintf('%s: %d of %d\n', mfilename, ii, num_steps);
end

% % <DEBUG>
% keyboard
% % </DEBUG>