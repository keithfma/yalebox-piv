function [rpk, cpk, stat] = yalebox_piv_peak_optim_lanczos(zz)
% function [rpk, cpk, stat] = yalebox_piv_peak_optim_lanczos(zz)
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

% init optimization
[rpk, cpk] = find(zz == max(zz(:)));
x0 = [rpk, cpk];

options = optimset('Display', 'none');

objective = @(x) -resample_lanczos(zz, x(1), x(2), 10);

% optimize
x = fmincon(objective, x0, [], [], [], [], x0-0.5, x0+0.5, [], options);

% return 
rpk = x(1);
cpk = x(2);
stat = 1;

end

function val = resample_lanczos(data, row, col, nkern)
% function val = resample_lanczos(data, row, col, nkern)
%
% Resample the input matrix using a lanczos (windowed sinc) interpolation
% kernel.
%
% Arguments:
%
%
%   nkern = Scalar, integer, half-width of interpolation kernel
%
% %

% extract data within kernel, ADD: with symmetric pad if needed
rr = (floor(row)-nkern+1):(floor(row)+nkern);
cc = (floor(col)-nkern+1):(floor(col)+nkern);
sub = data(rr, cc);

% compute kernel weights
cc = col-cc';
rr = row-rr;
kern = (sinc(cc).*sinc(cc/nkern))*(sinc(rr).*sinc(rr/nkern));
% kern = kern./sum(kern(:));

% resample
val = sum(sum(sub.*kern));

end