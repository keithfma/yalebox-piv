function [rpk, cpk, stat] = yalebox_piv_peak_optim_lanczos(zz, tol, maxitr)
% function [rpk, cpk, stat] = yalebox_piv_peak_optim_lanczos(zz, tol, maxitr)
%
% Find peak with subpixel accuracy using lanczos-kernel resampling and 2D
% optimization.
%
% Arguments:
%
% %

% set defaults
if nargin < 2; tol = 1e-6; end
if nargin < 3; maxitr = 100; end

% initialize optimization

% optimize


keyboard



end

function val = resample_lanczos(orig, row, col, nkern)
% function val = resample_lanczos(row, col, nkern)
%
% Resample the input matrix using a lanczos (windowed sinc) interpolation
% kernel.
%
% Arguments:
%
% %

% extract data within kernel, with symmetric pad if needed
r0 = floor(row);
c0 = floor(col);
rr = 

% compute kernel weights

% resample
val = []

end