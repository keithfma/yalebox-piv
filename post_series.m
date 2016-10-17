function [] = post_series(piv_netcdf, post_netcdf)
% 
% Run post-processing analyses on PIV data, and saves the results and
% metadata in a netCDF file.
%
% Arguments:
%   piv_netcdf: string, file name of netCDF containing PIV results as
%       produced by piv_series.m
%   post_netcdf: string, file name of netCDF to create for post-processing
%       results data and metadata
%
% Output netCDF is self-documenting
%
% % Keith Ma

% parse input arguments

% read in PIV data from netCDF

% create output file

% run analyses for each timestep

% finalize