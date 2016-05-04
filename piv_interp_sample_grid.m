function [r1, c1, u1, v1, x1, y1] = ...
    piv_interp_sample_grid(r0, c0, u0, v0, roi, len, spc, xw, yw, tension, verbose)
% function [r1, c1, u1, v1, x1, y1] = ...
%     piv_interp_sample_grid(r0, c0, u0, v0, roi, len, spc, xw, yw, tension, verbose)
%
% Get new sample grid and interpolate data to this new grid. Interpolation fills
% the entire new grid. For additinal piv passes, the new roi is imposed later by
% adding NaNs to the current estimate.
%
% Arguments:
%
%   r0, c0 = 2D matrix, original sample grid in pixel coordinates
%
%   u0, v0 = 2D matrix, displacement estimates on original sample grid
%
%   roi = 2D matrix, logical, flag indicating if point is within the data ROI
%       (1) or not (0)
%
%   len = Scalar, integer, side length for sample windows 
%
%   spc = Scalar, integer, spacing between sample windows
%   
%   xw, yw = Vectors, world coordinates at full image resolution
%
%   tension = Scalar, tension parameter for spline interpolation.
%
%   verbose = Logical flag, display verbose messages (1) or don't
% %
    
[r1, c1, x1, y1] = piv_sample_grid(len, spc, xw, yw);

if verbose
    fprintf('%s: new grid = [%d, %d], old grid = [%d, %d]\n', ...
        mfilename, size(r1,1), size(r1,2), size(r0,1), size(r0,2));
end

u1 = spline2d(c1(:), r1(:), c0(roi), r0(roi), u0(roi), tension);
v1 = spline2d(c1(:), r1(:), c0(roi), r0(roi), v0(roi), tension);

u1 = reshape(u1, size(r1));
v1 = reshape(v1, size(r1));