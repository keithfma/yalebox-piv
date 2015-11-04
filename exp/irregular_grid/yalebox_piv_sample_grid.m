function [rvec, cvec] = yalebox_piv_sample_grid(len, spc, sz)
% Create sample grid such that sample window edges fall at integer positions,
% and the grid is centered. Note that the first sample window requires 'len'
% rows and cols, and each additional window requires 'spc' rows and cols.
%
% EXPERIMENT: add an extra window so that the bottom sample window extends
%   beyond the image (and thus includes some padding)
%
% Arguments:
%
%   len = Scalar, integer, length of sample window in pixels
%
%   spc = Scalar, integer, grid spacing in pixels
%
%   sz = Vector, length == 2, grid dimensions [rows, columns] for the
%       original input data matrices (e.g. ini and fin).
%
%   rvec, cvec = Vector, integer, coordinate vectors for the sample grid in the
%       y (a.k.a. row) and x (a.k.a. column) directions, in pixels
% %

rem = mod(sz-len, spc)/2;

samp0 = floor(rem)+(len-1)/2; 
samp1 = sz-(ceil(rem)+(len-1)/2); 

% experiment {
samp0 = samp0-spc;
samp1 = samp1+spc;
% } experiment

rvec = samp0(1):spc:samp1(1);
cvec = samp0(2):spc:samp1(2);