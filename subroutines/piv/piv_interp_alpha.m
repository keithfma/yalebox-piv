function [u_to, v_to, roi_to] = piv_interp_alpha(...
    x_from, y_from, u_from, v_from, x_to, y_to, alpha, interp_method, ...
    extrap_method)
%
% function [u_to, v_to, roi_to] = piv_interp_alpha(...
%     x_from, y_from, u_from, v_from, x_to, y_to, alpha, interp_method, ...
%     extrap_method)
%
% Interpolate scattered points to regular grid using triangulation-based
% interpolant. Avoid known problems with "skinny" triangles along the data
% convex hull by applying different interpolation methods inside and
% outside of the data *alpha* hull. Alpha shapes are useful here because
% they are a subset of the Delaunay triangulation that excludes suspect
% triangles at the margins. Note that the interpolant (MATLAB
% scatteredInterpolant) is also based on a Delaunay triangulation.
%
% Arguments:
%   x_from, y_from = Location of scattered input points
%   u_from, v_from = Components of displacement vectors at scattered input points,
%       NaN values are ignored
%   x_to, y_to = Location of output points
%   alpha = Scalar, alpha shape radius parameter, in pixels, reccommended
%       to use 5*(sample window spacing), intrinsic (pixel) units
%   interp_method: Interpolation method applied inside the data alpha hull,
%       see scatteredInterpolant help for options, default is 'natural' if
%       empty or not set
%   extrap_method: Extrapolation method applied outside the data alpha
%       hull, see scatteredInterpolant help for options, default is
%       'nearest' if empty or not set
%
%   u_to, v_to = Interpolated output vectors where roi_to==1, NaNs elsewhere
%   roi_to = Region-of-interest mask for output, true elements are output. If a
%       scalar true is provided, the whole domain is used.
%
% References:
%   TODO: add alpha shape references (Edelsbrunner, etc)
% %

% parse input arguments
if nargin < 8 || isempty(interp_method)
    interp_method = 'natural';
end
if nargin < 9 || isempty(extrap_method)
    extrap_method = 'nearest';
end
    
% allocate outputs
u_to = nan(size(x_to));
v_to = nan(size(x_to));

% estimate ROI for interpolation from alpha hull of input points
shp = alphaShape(x_from(:), y_from(:), alpha);
roi_to = inShape(shp, x_to, y_to);

% interpolate within ROI
% NOTE: convert displacement vector to single complex number
% NOTE: extrapolation disabled - ROI points are all in convex hull
si = scatteredInterpolant(...
    x_from(:), y_from(:), u_from + 1i*v_from, interp_method);
tmp = si(x_to(roi_to), y_to(roi_to));
u_to(roi_to) = real(tmp);
v_to(roi_to) = imag(tmp);

% extrapolate outside ROI
si.Method = extrap_method;
si.ExtrapolationMethod = extrap_method;
tmp = si(x_to(~roi_to), y_to(~roi_to));
u_to(~roi_to) = real(tmp);
v_to(~roi_to) = imag(tmp);
