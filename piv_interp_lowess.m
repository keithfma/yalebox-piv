function [ug, vg] = piv_interp_lowess(xp, yp, up, vp, xg, yg, roi, npts, verbose)
% function [ug, vg] = piv_interp_lowess(xp, yp, up, vp, xg, yg, roi, npts, verbose)
%
% Smooth and interpolate scattered vectors to a regular grid using robust
% LOWESS. NaNs in input vector grids are ignored. Output vector grids are
% populated in the region-of-interest (roi) and NaN elsewhere.
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
%   npts = Scalar, number of points to include in local fit
%
%   ug, vg = 2D matrices, interpolated output vectors where roi==1, NaNs elsewhere
%
%   verbose = Scalar, logical, display verbose messages (1) or don't (0)
% %

% constants
robust_opt = 'bisquare';

% init
from = ~isnan(up) & ~isnan(vp);
span = npts/sum(from(:));

if verbose
    fprintf('%s: span = %d pts (%.3f%%), robust = %s, source pts = %d, dest pts = %d\n', ...
        mfilename, npts, span, robust_opt, sum(from(:)), sum(roi(:)));
end

% interpolate
ug = nan(size(roi));
u_model = fit([xp(from), yp(from)], up(from), 'lowess', 'Span', span, 'Robust', robust_opt);
ug(roi) = u_model(xg(roi), yg(roi));

vg = nan(size(roi));
v_model = fit([xp(from), yp(from)], vp(from), 'lowess', 'Span', span, 'Robust', robust_opt);
vg(roi) = v_model(xg(roi), yg(roi));
