function [r_peak, c_peak, val_peak] = piv_peak_optim_interp(zz, precision)
%
% NEW: find peak by optimization of interpolated surface, terminate at precision
%
% Locate subpixel peak with specified precision by upscaling in the vicinity of
% the peak. Uses MATLAB's built in interpolation routines for simplicity.
%
% Arguments:
%
%   zz = 2D matrix, data set for which to estimate peak
%
%   precision = Scalar, grid spacing in fine (high-res) grid, which gives the
%       precision of the subplixel peak location estimate
% %

% check for sane inputs
validateattributes(zz, {'numeric'}, {'2d', 'real'});
validateattributes(precision, {'numeric'}, {'scalar', '<', 1});

% initial guess at maximum
% TODO: use max to get index, then convert to subscript
[r_peak_init, c_peak_init] = find(zz == max(zz(:)));

% fail if no unique maximum 
if numel(r_peak_init) ~= 1  
    warning('Failed to find subpixel peak');
    r_peak = NaN;
    c_peak = NaN;
    val_peak = NaN;
    return
end

% define the objective function
[rr, cc] = ndgrid(1:size(zz, 1), 1:size(zz, 2));
gi = griddedInterpolant(rr, cc, -zz, 'spline'); % negate zz to so peak -> minimum
objective = @(x) gi(x(1), x(2));

% find peak by optimization 
% NOTE: The options are probably over-specified
opt = optimoptions(@fminunc);
opt.Algorithm = 'quasi-newton';
opt.Display = 'notify';
opt.FunctionTolerance = 0; % terminate based on peak change, not value change
opt.StepTolerance = precision;
opt.FiniteDifferenceType = 'central';
opt.TypicalX = [r_peak_init, c_peak_init];
[tmp, val_peak] = fminunc(objective, [r_peak_init, c_peak_init], opt);
r_peak = tmp(1);
c_peak = tmp(2);

end
