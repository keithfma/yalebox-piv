function [rpk, cpk, val] = piv_peak_interp(zz, precision)
%
% Locate subpixel peak with specified precision by brute-force interpolation.
% %

% debug: hard-coded inputs for testing
zz = peaks(100);
precision = 0.001;

% assume sane inputs for performance

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



function [rr, cc, vv, ss] = fail()
% Set output values and throw a warning if the main function fails to find the
% subpixel peak.

    warning('Failed to find subpixel peak');
    rr = NaN;
    cc = NaN;
    vv = NaN;
    
end