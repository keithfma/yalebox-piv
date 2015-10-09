function [rpk, cpk, stat] = yalebox_piv_peak_optim_lanczos(zz, tol, maxitr)
% function [rpk, cpk, stat] = yalebox_piv_peak_optim_lanczos(zz, tol, maxitr)
%
% Find peak with subpixel accuracy using lanczos-kernel resampling and 2D
% optimization.

% set defaults
if nargin < 2; tol = 1e-6; end
if nargin < 3; maxitr = 100; end

keyboard



