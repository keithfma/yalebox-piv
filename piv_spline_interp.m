function [ug, vg] = piv_spline_interp(xp, yp, up, vp, xg, yg, roi, tension, verbose)
% function [ug, vg] = piv_spline_interp(xp, yp, up, vp, xg, yg, roi, tension, verbose)
%
% Smooth and interpolate scattered vectors to a regular grid using
% spline-in-tension interpolation (see [1]). NaNs in input grids are ignored. Output grids
% are populated in the region-of-interest (roi) and NaN elsewhere.
%
% Arguments:
%
%   xp, yp = 2D matrices, location of scattered input points
%
%   up, vp = 2D matrices, components of displacement vectors at scattered input 
%       points, NaN values are ignored
%
%   xg, yg = 2D matrices, regular grid for output vectors
%
%   roi = 2D matrix, region-of-interest mask, 1 vectors should be output
%
%   tension = Scalar, tension parameter for interpolation routine, see [1]
%
%   ug, vg = 2D matrices, interpolated output vectors where roi==1, NaNs elsewhere
%
%   verbose = Scalar, logical, display verbose messages (1) or don't (0)
%
% References:
% [1]
% 
% %

% init
from = ~isnan(up) & ~isnan(vp);

if verbose
    fprintf('%s: source pts = %d, dest pts = %d, tension = %f\n',...
        mfilename, sum(from(:)), sum(roi(:)), tension);
end

% interpolate
ug = nan(size(roi));
ug(roi) = spline2d(xg(roi), yg(roi), xp(from), yp(from), up(from), tension);

vg = nan(size(roi));
vg(roi) = spline2d(xg(roi), yg(roi), xp(from), yp(from), vp(from), tension);