function [u_out, v_out] = piv_interp_linear_alpha(...
    x_in, y_in, u_in, v_in, x_out, y_out, roi_out, alpha, verbose)
%
% TODO: add an ROI for input data
%
% Interpolate scattered vectors using linear interpolation and extrapolation
% based on Delaunay triangulation (MATLAB scatteredInterpolant). NaNs in input
% vector grids are ignored. Output vector grids are populated in the
% region-of-interest (roi) *and* within the alpha hull of the input data, and
% NaN elsewhere.
%
% *NEW* Points outside the alpha shape are extrapolated, not interpolated. This
% avoids the nefarious influence of skinny triangles connecting the concave
% boundaries of the wedge.
%
% Since the alpha shape is closely related to the Delaunay triangulation...
%
% TODO: update defs
% Arguments:
%   xi, yi = Location of scattered input points
%   uin, vin = Components of displacement vectors at scattered input points,
%       NaN values are ignored
%   xout, yout = Location of output points
%   roi = Region-of-interest mask for output, true elements are output. If a
%       scalar true is provided, the whole domain is used.
%   alpha = Scalar, alpha shape radius parameter 
%   uout, vout = Interpolated output vectors where roi==1, NaNs elsewhere
%   verbose = Scalar, logical, display verbose messages (1) or don't (0)
% % 

% init
from = ~isnan(u_in) & ~isnan(v_in);
if isscalar(roi_out)
    if roi_out
        roi_out = true(size(x_out));
    else
        error('%s: Bad value for argument "roi_out"', mfilename);
    end
end

% create interpolation / extrapolation masks using alpha shape
% trying alpha = 3*sampspc...
alpha_shp = alphaShape(x_in(from), y_in(from), alpha);
to = roi_out & inShape(alpha_shp, x_out, y_out);

if verbose
    fprintf('%s: source pts = %d, dest pts = %d\n', ...
        mfilename, sum(from(:)), sum(to(:)));
end

% interpolate vectors as complex field
c_out = zeros(size(to));
interpolant = scatteredInterpolant(x_in(from), y_in(from), ...
    u_in(from) + 1i*v_in(from), 'linear', 'linear');
c_out(to) = interpolant(x_out(to), y_out(to));

% unpack complex values and re-apply mask
% EXPERIMENT: set out-of-roi values to 0, not NaN
u_out = real(c_out);
u_out(~to) = 0;

v_out = imag(c_out);
v_out(~to) = 0;