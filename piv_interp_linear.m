function [ug, vg] = piv_interp_linear(xp, yp, up, vp, xg, yg, roi, verbose)
% function [ug, vg] = piv_interp_linear(xp, yp, up, vp, xg, yg, roi, verbose)
%
% Interpolate scattered vectors to a regular grid using linear interpolation and
% extrapolation based on Delaunay triangulation (MATLAB scatteredInterpolant).
% NaNs in input vector grids are ignored. Output vector grids are populated in
% the region-of-interest (roi) and NaN elsewhere.
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
%   ug, vg = 2D matrices, interpolated output vectors where roi==1, NaNs elsewhere
%
%   verbose = Scalar, logical, display verbose messages (1) or don't (0)
% %

% init
from = ~isnan(up) & ~isnan(vp);

if verbose
    fprintf('%s: source pts = %d, dest pts = %d\n', ...
        mfilename, sum(from(:)), sum(roi(:)));
end

% interpolate
ug = nan(size(roi));
interpolant = scatteredInterpolant(xp(from), yp(from), up(from), 'linear', 'linear');
ug(roi) = interpolant(xg(roi), yg(roi));

vg = nan(size(roi));
interpolant.Values = vp(from);
vg(roi) = interpolant(xg(roi), yg(roi));
