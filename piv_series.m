function [] = piv_series(output_file, input_file, samplen, sampspc, intrlen, ...
                  npass, valid_max, valid_eps, lowess_span_pts, spline_tension, ...
                  min_frac_data, min_frac_overlap, verbose)
% function [] = piv_series(output_file, input_file, samplen, sampspc, intrlen, ...
%                   npass, valid_max, valid_eps, lowess_span_pts, spline_tension, ...
%                   min_frac_data, min_frac_overlap, verbose)
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
% verbose = Scalar, logical, display verbose messages for this function and its
%   children (1) or don't (0) 
% %

% check for sane arguments (pass-through arguments are checked in subroutines)
narginchk(13, 13); 
validateattributes(output_file, {'char'}, {'vector'});
validateattributes(input_file, {'char'}, {'vector'});

% read input dimension values
x_img =      double( ncread(input_file, 'x')     );
y_img =      double( ncread(input_file, 'y')     );
step_img =   double( ncread(input_file, 'step')  );

% compute output dimensions
[~, ~, x_piv, y_piv] = piv_sample_grid(samplen(end), sampspc(end), x_img, y_img);        
step_piv = step_img(1:end-1)+0.5;
nx = length(x_piv);
ny = length(y_piv);
ns = length(step_piv);

% create netcdf file
ncid = netcdf.create(output_file, 'NETCDF4');

% add global attributes
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv samplen', samplen);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv sampspc', sampspc);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv intrlen', intrlen);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv npass', npass);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv valid_max', valid_max);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv valid_eps', valid_eps);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv valid_eps', valid_eps);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv lowess_span_pts', lowess_span_pts);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv spline_tension', spline_tension);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv min_frac_data', min_frac_data);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'piv min_frac_overlap', min_frac_overlap);
netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'git commit hash', util_git_hash());
%netcdf.putAtt(ncid, netcdf.getConstant('GLOBAL'), 'input file MD5 hash', util_md5_hash(input_file));

% create dimensions
x_dimid = netcdf.defDim(ncid, 'x', nx);
y_dimid = netcdf.defDim(ncid, 'y', ny);
s_dimid = netcdf.defDim(ncid, 'step', ns);

% define variables and thier attributes, compression, and chunking
dim_3d = [y_dimid, x_dimid, s_dimid];
chunk_3d = [length(y_piv), length(x_piv), 1];

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
netcdf.putVar(ncid, x_varid, x_piv);
netcdf.putVar(ncid, y_varid, y_piv);
netcdf.putVar(ncid, s_varid, step_piv);
netcdf.close(ncid);

% initialize loop by reading first image
roi_const = logical(ncread(input_file, 'mask_manual', [1, 1], [inf, inf]));
roi1 = logical(ncread(input_file, 'mask_auto', [1, 1, 1], [inf, inf, 1]))& roi_const;
% img1 = double(ncread(input_file, 'img', [1, 1, 1], [inf, inf, 1]));
img1 = double(ncread(input_file, 'img_rgb', [1, 1, 1, 1], [inf, inf, inf, 1]));
img1 = rgb2gray(img1./255);
img1(~roi1) = 0;

% analyse all steps
for ii = 1:ns
    
    if verbose
        fprintf('\n%s: begin step = %.1f\n', mfilename, step_piv(ii));
    end
    
    % update image and roi pair
    img0 = img1;
    roi0 = roi1;
    
    roi1 = logical(ncread(input_file, 'mask_auto', [1, 1, ii+1], [inf, inf, 1])) & roi_const;
    % img1 = double(ncread(input_file, 'img', [1, 1, ii+1], [inf, inf, 1]));
    img1 = double(ncread(input_file, 'img_rgb', [1, 1, 1, ii+1], [inf, inf, inf, 1]));
    img1 = rgb2gray(img1./255);
    img1(~roi1) = 0;
    
    % perform piv analysis
    [~, ~, u_piv, v_piv, roi_piv] = ...
        piv(img0, img1, roi0, roi1, x_img, y_img, samplen, sampspc, intrlen, npass, ...
            valid_max, valid_eps, lowess_span_pts, spline_tension, ...
            min_frac_data, min_frac_overlap, verbose); 
        
    % write results to output file
    ncid = netcdf.open(output_file, 'WRITE');    
    netcdf.putVar(ncid, u_varid, [0, 0, ii-1], [ny, nx, 1], u_piv);
    netcdf.putVar(ncid, v_varid, [0, 0, ii-1], [ny, nx, 1], v_piv); 
    netcdf.putVar(ncid, r_varid, [0, 0, ii-1], [ny, nx, 1], int8(roi_piv));  
    netcdf.close(ncid);
    
    if verbose
        fprintf('%s: end step = %.1f\n', mfilename,  step_piv(ii));
    end
    
end
