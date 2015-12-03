function [r_cntr, c_cntr] = yalebox_piv_centroid(r0, c0, window, mask_val)
% function [r_cntr, c_cntr] = yalebox_piv_centroid(r0, c0, window, mask_val)
%
% Find and return the centroid of the data (i.e. sand) in the input window.
%
% Arguments:
%
% r0, c0 = Scalar, the row, column index in the parent matrix that corresponds
%   to element (1,1) in 'window'
%
% window = 2D matrix, input data for which to find the centroid.
%
% mask_val = Scalar, value assigned to "no data" positions in 'window'
% %

% compute centroid
[r_data, c_data] = find(window ~= mask_val);
n_data = length(r_data);
r_cntr = sum(r_data)/n_data;
c_cntr = sum(c_data)/n_data;

% convert to location in parent matrix
r_cntr = r_cntr+r0-1;
c_cntr = c_cntr+c0-1;