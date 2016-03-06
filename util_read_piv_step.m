function [step, xx, yy, uu, vv, mm] = util_read_piv_step(piv_file, index)
%
% Reads single step from a PIV results netCDF file as produced by piv_series()
%
% Arguments:
%
% piv_file = String, file name of netCDF file containing output from the
%   piv_series() analysis
%
% index = Scalar, index of timestep in piv_file to plot, uses MATLAB-style
%   1-based indices.
%
% step = 
%
% xx, yy = 
%
% uu, vv, mm = 
% %

% read data from netCDF
xx = ncread(piv_file, 'x'); 
xx = double(xx);
nx = numel(xx);

yy = ncread(piv_file, 'y'); 
yy = double(yy);
ny = numel(yy);

step = ncread(piv_file, 'step', index, 1);
step = double(step);

uu = ncread(piv_file, 'u', [1, 1, index], [nx, ny, 1])';
uu = double(squeeze(uu));

vv = ncread(piv_file, 'v', [1, 1, index], [nx, ny, 1])';
vv = double(squeeze(vv));

mm = sqrt(uu.*uu+vv.*vv);