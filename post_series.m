function [] = post_series(piv_netcdf, post_netcdf)
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
%
% Output netCDF is self-documenting
%
% % Keith Ma

% parse input arguments
validateattributes(piv_netcdf, {'char'}, {'nonempty'}, mfilename, 'piv_netcdf');
assert(exist(piv_netcdf, 'file') == 2, ...
    sprintf('input file %s does not exist', piv_netcdf));

validateattributes(post_netcdf, {'char'}, {'nonempty'}, mfilename, 'post_netcdf');
assert(exist(post_netcdf, 'file') ~= 2, ...
    sprintf('output file %s already exists', piv_netcdf));    

% read in PIV data from netCDF
xx   = double(ncread(piv_netcdf, 'x')); 
yy   = double(ncread(piv_netcdf, 'y')); 
step = double(ncread(piv_netcdf, 'step'));
uu   = double(ncread(piv_netcdf, 'u'));
vv   = double(ncread(piv_netcdf, 'v'));
roi  = double(ncread(piv_netcdf, 'roi'));

% create output file

% run analyses for each timestep
% % pro- and retro-plate velocity
% % area budget
% % infinitesimal strain analysis
% % finite strain analysis

% finalize

% debug
keyboard