function [uout, vout] = piv_interp_linear(...
    xin, yin, uin, vin, xout, yout, roi, verbose)
% function [uout, vout] = piv_interp_linear(...
%   xin, yin, uin, vin, xout, yout, roi, verbose)
%
% Interpolate scattered vectors using linear interpolation and extrapolation
% based on Delaunay triangulation (MATLAB scatteredInterpolant). NaNs in input
% vector grids are ignored. Output vector grids are populated in the
% region-of-interest (roi) and NaN elsewhere.
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
    fprintf('%s: source pts = %d, dest pts = %d\n', ...
        mfilename, sum(from(:)), sum(roi(:)));
end

% interpolate vectors as complex field
cout = nan(size(roi));
interpolant = scatteredInterpolant(xin(from), yin(from), ...
    uin(from) + 1i*vin(from), 'linear', 'linear');
cout(roi) = interpolant(xout(roi), yout(roi));
uout = real(cout);
vout = imag(cout);