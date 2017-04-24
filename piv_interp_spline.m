function [uout, vout] = piv_interp_spline(...
    xi, yi, uin, vin, xout, yout, roi, tension, verbose)
% function [uout, vout] = piv_interp_spline(...
%     xi, yi, uin, vin, xout, yout, roi, tension, verbose)
%
% Smooth and interpolate scattered vectors within ROI using spline-in-tension
% interpolation (see [1]). NaNs in input grids are ignored. Output grids are
% populated in the region-of-interest (roi) and NaN elsewhere.
%
% Arguments:
%   xi, yi = Location of scattered input points
%   uin, vin = Components of displacement vectors at scattered input points,
%       NaN values are ignored
%   xout, yout = Location of output points
%   roi = Region-of-interest mask for output, true elements are output. If a
%       scalar true is provided, the whole domain is used.
%   tension = Scalar, tension parameter for interpolation routine, see [1]
%   uout, vout = Interpolated output vectors where roi==1, NaNs elsewhere
%   verbose = Scalar, logical, display verbose messages (1) or don't (0)
%
% References: 
% [1] Wessel, P., & Bercovici, D. (1998). Interpolation with splines in tension:
%     A Green’s function approach. Mathematical Geology, 30(1), 77–93.
% %

% init
from = ~isnan(uin) & ~isnan(vin);
if isscalar(roi)
    if roi
        roi = true(size(xout));
    else
        error('%s: Bad value for argument "roi"', mfilename);
    end
end

if verbose
    fprintf('%s: source pts = %d, dest pts = %d, tension = %f\n',...
        mfilename, sum(from(:)), sum(roi(:)), tension);
end

% NOTE: ~2x faster to interpolate as a complex field
cout = nan(size(roi));
cout(roi) = spline2d(xout(roi), yout(roi), xi(from), yi(from), ...
                     uin(from) + 1i*vin(from), tension);
uout = real(cout);
vout = imag(cout);