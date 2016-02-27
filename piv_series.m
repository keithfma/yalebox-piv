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

% check for sane input file (only checks for names of dims and variables)
info = ncinfo(input_file);

% ...3 dimensions: x, y, step
assert( length(info.Dimensions) == 3 );
for ii = 1:3
    assert( ismember(info.Dimensions(ii).Name, {'x', 'y', 'step'}) );
end

% ...6 variables: x, y, step, intensity, mask_manual, mask_auto
assert( length(info.Variables) == 6 );
for ii = 1:6
    assert( ismember(info.Variables(ii).Name, {'x', 'y', 'step', 'intensity', ...
        'mask_manual', 'mask_auto'}) );
end

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

% define variables and thier attributes, compression, and chunking

% finish netcdf creation
netcdf.endDef(ncid);
netcdf.close(ncid);

% populate constant variables

% loop over all timesteps

    % perform piv analysis

% done