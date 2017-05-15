function [u_to, v_to] = piv_interp_spline(...
    x_from, y_from, u_from, v_from, x_to, y_to, roi_to, tension, verbose)
% function [u_to, v_to] = piv_interp_spline(...
%     x_from, y_from, u_from, v_from, x_to, y_to, roi_to, tension, verbose)
%
% Smooth and interpolate scattered vectors within ROI using spline-in-tension
% interpolation (see [1]). NaNs in input grids are ignored. Output grids are
% populated in the region-of-interest (roi_to) and NaN elsewhere.
%
% Arguments:
%   x_from, y_from = Location of scattered input points
%   u_from, v_from = Components of displacement vectors at scattered input points,
%       NaN values are ignored
%   x_to, y_to = Location of output points
%   roi_to = Region-of-interest mask for output, true elements are output. If a
%       scalar true is provided, the whole domain is used.
%   tension = Scalar, tension parameter for interpolation routine, see [1]
%   u_to, v_to = Interpolated output vectors where roi_to==1, NaNs elsewhere
%   verbose = Scalar, logical, display verbose messages (1) or don't (0)
%
% References: 
% [1] Wessel, P., & Bercovici, D. (1998). Interpolation with splines in tension:
%     A Green?s function approach. Mathematical Geology, 30(1), 77?93.
% %

% parse roi input arguments
from = ~isnan(u_from);
if isscalar(roi_to) && roi_to
    roi_to = true(size(x_to));
end
    
if verbose
    fprintf('%s: source pts = %d, dest pts = %d, tension = %f\n',...
        mfilename, sum(from(:)), sum(roi_to(:)), tension);
end

% interpolate as a complex field (~2x faster)
c_to = zeros(size(roi_to));
c_to(roi_to) = spline2d(x_to(roi_to), y_to(roi_to), ...
                        x_from(from), y_from(from), ...
                        u_from(from) + 1i*v_from(from), ...
                        tension);
u_to = real(c_to);
v_to = imag(c_to);

% re-apply ROI (necessary because no such thing as imaginary NaN)
u_to(~roi_to) = NaN;
v_to(~roi_to) = NaN;
