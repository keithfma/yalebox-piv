function [rpk, cpk, stat] = yalebox_piv_peak_radsym(zz)
% function [rpk, cpk, stat] = yalebox_piv_peak_radsym(zz)
%
% Find the position of the peak in matrix zz with subpixel accuracy. Peak
% location is determined from an explicit solution the best-fit center of
% radial symmetry. See reference [1] for details and the core function
% radialcenter(). Solution is computed on a square subset of the data
% around the integer-precision peak. This window size is hard-coded at
% present.
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
%   stat = Logical, scalar, return status flag, true if successful, false if
%       unsuccessful
%
% References:
%
% [1] Parthasarathy, R. (2012). Rapid, accurate particle tracking by
% calculation of radial symmetry centers. Nature Methods, 9(7), 724?726.
% doi:10.1038/nmeth.2071

% define parameters
nwin = 7; % odd integer, square dimension of the subset window used to find the peak
delta = (nwin-1)/2; 

% find integer peak position
[rpk, cpk] = find(zz == max(zz(:)));

% check failure conditions
%   1-2) no unique maximum
%   3-6) peak window extends outside the edge of the matrix
stat = true;
if numel(rpk) ~= 1 || ...
        numel(cpk) ~= 1 || ...
        rpk < 1+delta || ...
        rpk > size(zz, 1)-delta || ...
        cpk < 1+delta || ...
        cpk > size(zz, 2)-delta
    stat = false;
    return
end
    
% extract peak window
peak_win = zz((rpk-delta):(rpk+delta), (cpk-delta):(cpk+delta));
                     
% compute sub-pixel peak location in window
[cpk_win, rpk_win] = radialcenter(peak_win);

% convert to peak location in zz
cpk = cpk_win-delta-1+cpk;
rpk = rpk_win-delta-1+rpk;
