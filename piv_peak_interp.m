function [rpk, cpk, val] = piv_peak_interp(zz, precision)
%
% Locate subpixel peak with specified precision by upscaling in the vicinity of
% the peak. Uses MATLAB's built in interpolation routines for simplicity.
%
% Arguments:
%
%   zz = 2D matrix, data set for which to estimate peak
%
%   precision = Scalar, grid spacing in fine (high-res) grid, which gives the
%       precision of the subplixel peak location estimate, must divide evenly
%       into 1 (e.g 0.1 is OK, 0.3 is not) so that the fine grid is aligned with
%       the coarse.
% %

% check for sane inputs
validateattributes(zz, {'numeric'}, {'2d', 'real'});
validateattributes(precision, {'numeric'}, {'scalar', '<', 1});
assert(mod(1/precision, 1) == 0);

% get peak, fail if no unique maximum 
[rpk, cpk] = find(zz == max(zz(:)));
if numel(rpk) ~= 1  
    [rpk, cpk, val] = fail();
    return
end

% generate grids: coarse for zz, fine in the vicnity of the peak w/ no extrapolation
[z_num_row, z_num_col] = size(zz);
row_fine0 = max(rpk-1, 1);
col_fine0 = max(cpk-1, 1);
row_fine1 = min(rpk+1, z_num_row);
col_fine1 = min(cpk+1, z_num_col);
[col, row] = meshgrid(1:z_num_col, 1:z_num_row);
[col_fine, row_fine] = meshgrid(col_fine0:precision:col_fine1, row_fine0:precision:row_fine1); 

% interpolate and extract peak location and value
zz_fine = interp2(col, row, zz, col_fine, row_fine, 'spline');
val = max(zz_fine(:));
ipk = find(zz_fine == val);
cpk = col_fine(ipk);
rpk = row_fine(ipk);

end

function [rr, cc, vv] = fail()
% Set output values and throw a warning if the main function fails to find the
% subpixel peak.

    warning('Failed to find subpixel peak');
    rr = NaN;
    cc = NaN;
    vv = NaN;
    
end