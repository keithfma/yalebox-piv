function [] = piv_series(output_file, input_file, samplen, sampspc, intrlen, ...
                  npass, valid_max, valid_eps)
% 
% Run PIV analysis for a given input series. Input is expected to be a netCDF
% file as created by prep_series(). Results are saved in a new netCDF file which
% includes all relevant metadata.
%
% output_file = String, name of the output netCDF file containing PIV results
%
% input_file = String, name of the input netCDF file containing pre-processed
%   image data
%
% samplen, sampspc, intrlen, npass, valid_max, valid_eps = Select input
%   variables for PIV routine piv(), other inputs are contained in the input
%   file.
%
% %

% check for sane arguments (pass-through arguments are checked in subroutines)
% narginchk(8, 8); 
validateattributes(output_file, {'char'}, {'vector'});
validateattributes(input_file, {'char'}, {'vector'});

% check for sane input file
info = ncinfo(input_file);
assert( length(info.Dimensions) == 3 );
for ii = 1:3
    assert( ismember(info.Dimensions(ii).Name, {'x', 'y', 'step'}), ...
        sprintf('Invalid input file, failed to find dimension %s', info.Dimensions(ii).Name));
end
assert( length(info.Variables) == 6 );
for ii = 1:6
    assert( ismember(info.Variables(ii).Name, {'x', 'y', 'step', 'intensity', 'mask_manual', 'mask_auto'}), ...
              sprintf('Invalid input file, failed to find variable %s', info.Variables(ii).Name));
end

% read dimension values
xx = ncread(input_file, 'x');
yy = ncread(input_file, 'y');
step = ncread(input_file, 'step');
step = double(step);

% create netcdf file
ncid = netcdf.create(output_file, 'NETCDF4');

% add global attributes
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv samplen', samplen);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv sampspc', sampspc);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv intrlen', intrlen);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv npass', npass);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv valid_max', valid_max);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv valid_eps', valid_eps);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'git commit hash', util_git_hash());
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'input file MD5 hash', util_md5_hash(input_file));
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv valid_eps', valid_eps);

% create dimensions
x_dimid = netcdf.defDim(ncid, 'x', length(xx));
y_dimid = netcdf.defDim(ncid, 'y', length(yy));
s_dimid = netcdf.defDim(ncid, 'step', length(step)-1);

% define variables and thier attributes, compression, and chunking
dim_3d = [x_dimid, y_dimid, s_dimid];
chunk_3d = [length(xx), length(yy), 1];

x_varid = netcdf.defVar(ncid, 'x', 'NC_FLOAT', x_dimid);
netcdf.putAtt(ncid, x_varid, 'long_name', 'horizontal position');
netcdf.putAtt(ncid, x_varid, 'units', 'meters');

y_varid = netcdf.defVar(ncid, 'y', 'NC_FLOAT', y_dimid);
netcdf.putAtt(ncid, y_varid, 'long_name', 'vertical position');
netcdf.putAtt(ncid, y_varid, 'units', 'meters');

s_varid = netcdf.defVar(ncid, 'step', 'NC_FLOAT', s_dimid);
netcdf.putAtt(ncid, s_varid, 'long_name', 'step number');
netcdf.putAtt(ncid, s_varid, 'units', '1');

u_varid = netcdf.defVar(ncid, 'u', 'NC_FLOAT', dim_3d);
netcdf.putAtt(ncid, u_varid, 'long_name', 'displacement vector, x-component');
netcdf.putAtt(ncid, u_varid, 'units', 'meters/step');
netcdf.defVarDeflate(ncid, u_varid, true, true, 1);
netcdf.defVarChunking(ncid, u_varid, 'CHUNKED', chunk_3d);

v_varid = netcdf.defVar(ncid, 'v', 'NC_FLOAT', dim_3d);
netcdf.putAtt(ncid, v_varid, 'long_name', 'displacement vector, y-component');
netcdf.putAtt(ncid, v_varid, 'units', 'meters/step');
netcdf.defVarDeflate(ncid, v_varid, true, true, 1);
netcdf.defVarChunking(ncid, v_varid, 'CHUNKED', chunk_3d);

r_varid = netcdf.defVar(ncid, 'roi', 'NC_BYTE', dim_3d);
netcdf.putAtt(ncid, r_varid, 'long_name', 'displacement vector mask');
netcdf.putAtt(ncid, r_varid, 'units', 'boolean');
netcdf.defVarDeflate(ncid, r_varid, true, true, 1);
netcdf.defVarChunking(ncid, r_varid, 'CHUNKED', chunk_3d);

% finish netcdf creation
netcdf.endDef(ncid);
netcdf.close(ncid);

% populate dimension values
ncid = netcdf.open(output_file, 'WRITE');
netcdf.putVar(ncid, x_varid, xx);
netcdf.putVar(ncid, y_varid, yy);
netcdf.putVar(ncid, s_varid, step(1:end-1)+0.5);
netcdf.close(ncid);

% analyse all steps
roi_const = ncread(input_file, 'mask_manual', [1, 1], [inf, inf])';
roi1 = ncread(input_file, 'mask_auto', [1, 1, 1], [inf, inf, 1])' & roi_const;
img1 = ncread(input_file, 'intensity', [1, 1, 1], [inf, inf, 1])';
for ii = 1:length(step)-1
    
    % update image and roi pair
    img0 = img1;
    roi0 = roi1;
    img1 = ncread(input_file, 'intensity', [1, 1, ii+1], [inf, inf, 1])';
    roi1 = ncread(input_file, 'mask_auto', [1, 1, ii+1], [inf, inf, 1])' & roi_const;
    
    % perform piv analysis
    [x_piv, y_piv, u_piv, v_piv, roi_piv] = ...
        piv(double(img0), double(img1), roi0, roi1, double(xx), double(yy), samplen, sampspc, intrlen, npass, ...
            valid_max, valid_eps, 1); 
        
    % debug
    subplot(2,1,1); imagesc(u_piv);
    subplot(2,1,2); imagesc(v_piv);
    pause;
    
    
    % write results to output file

end

%% subroutines

function img = fetch_img(file, index)
%
% Read image data from the input netCDF file for the specified time step


