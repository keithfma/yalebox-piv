function [ug, vg] = piv_lowess_interp(xp, yp, up, vp, xg, yg, roi, npts)
% function [ug, vg] = piv_lowess_interp(xp, yp, up, vp, xg, yg, roi, npts)
%
% Smooth and interpolate scattered vectors to a regular grid using robust
% LOWESS. NaNs in input vector grids are ignored. Output vector grids are
% populated in the region-of-interest (roi) and NaN elsewhere.
%
% Arguments:
%   xp, yp = 2D matrices, location of scattered input points
%   up, vp = 2D matrices, components of displacement vectors at scattered input 
%       points, NaN values are ignored
%   xg, yg = 2D matrices, regular grid for output vectors
%   roi = 2D matrix, region-of-interest mask, 1 vectors should be output
%   npts = Scalar, number of points to include in local fit
%   ug, vg = 2D matrices, interpolated output vectors where roi==1, NaNs elsewhere
% 
% %

from = ~isnan(up) & ~isnan(vp);
span = npts/sum(from(:));

ug = nan(size(roi));
u_model = fit([xp(from), yp(from)], up(from), 'lowess', 'Span', span, 'Robust', 'bisquare');
ug(roi) = u_model(xg(roi), yg(roi));

vg = nan(size(roi));
v_model = fit([xp(from), yp(from)], vp(from), 'lowess', 'Span', span, 'Robust', 'bisquare');
vg(roi) = v_model(xg(roi), yg(roi));

end