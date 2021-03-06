function [rpk, cpk, val] = piv_peak_gauss2d(zz)
% function [rpk, cpk, val] = piv_peak_gauss2d(zz)
%
% Find the position of the peak in matrix zz with subpixel accuracy. Peak
% location is determined from an explicit solution of two-dimensional
% Gaussian regression [1]. If the peak cannot be fit at subpixel
% accuracy, no peak is returned (see Arguments). This choice reflects the
% fact that a lack of subpixel displacement causes spurious gradients - it
% is preferable to drop the vector and interpolate.
%
% Arguments:
%   zz = 2D matrix, data plane in which to locate the peak
%
%   rpk = Scalar, double, row-coordinate location of the peak, set to -1 if
%       the peak cannot be fit.
%
%   cpk = Scalar, double, column-coordinate location of the peak, set to -1
%       if the peak cannot be fit
%
%   val = Scalar, double, peak value from the best-fit gaussian
%
% References:
%
% [1] Nobach, H., & Honkanen, M. (2005). Two-dimensional Gaussian regression for
% sub-pixel displacement estimation in particle image velocimetry or particle
% position estimation in particle tracking velocimetry. Experiments in Fluids,
% 38(4), 511–515. doi:10.1007/s00348-005-0942-3

[rpk, cpk] = find(zz == max(zz(:)));

% fail if no unique maximum
if numel(rpk) ~= 1 || numel(cpk) ~= 1 
    [rpk, cpk, val] = fail();
    return
end

% fail if peak is at the edge of the matrix
if rpk == 1 || rpk == size(zz, 1) || cpk == 1 || cpk == size(zz,2)
    [rpk, cpk, val] = fail();
    return
end

% offset to eliminate non-positive (gaussian is always positive)
offset = min(zz(:));
zz = zz-offset+eps;

% compute coefficients
c10 = 0;
c01 = 0;
c11 = 0;
c20 = 0;
c02 = 0;
c00 = 0;
for ii = -1:1
    for jj = -1:1
        logterm = log(zz(rpk+jj,cpk+ii));
        c10 = c10 + ii*logterm/6;
        c01 = c01 + jj*logterm/6;
        c11 = c11 + ii*jj*logterm/4;
        c20 = c20 + (3*ii^2-2)*logterm/6;
        c02 = c02 + (3*jj^2-2)*logterm/6;
        c00 = c00 + (5-3*ii^2-3*jj^2)*logterm/9;
    end
end

% compute sub-pixel displacement
dr = ( c11*c10-2*c01*c20 )/( 4*c20*c02 - c11^2 );
dc = ( c11*c01-2*c10*c02 )/( 4*c20*c02 - c11^2 );

% compute peak value
val = exp( c00-c20*dc^2-c11*dc*dr-c02*dr^2 )+offset;

% fail if subpixel fit erroneously large
if abs(dr) >= 1 || abs(dc) >= 1
    [rpk, cpk, val] = fail();
end

% apply subpixel displacement
rpk = rpk+dr;
cpk = cpk+dc;
    
end

function [rr, cc, vv, ss] = fail()
% Set output values and throw a warning if the main function fails to find the
% subpixel peak.

    % warning('Failed to find subpixel peak'); % NOTE: this turns out to be more annoying than helpful
    rr = NaN;
    cc = NaN;
    vv = NaN;
    
end
