function [rpk, cpk, stat] = yalebox_piv_peak_optim_sinc(cspace)
% function [rpk, cpk, stat] = yalebox_piv_peak_optim_sinc(cspace)
%
% Find peak with subpixel accuracy using sinc resampling in the fourier domain
% and 2D optimization.
%
% Arguments:
%
% %

% init optimization
[rpk, cpk] = find(cspace == max(cspace(:)));
x0 = [rpk, cpk];
cwave = fft2(cspace);
keyboard

% options = optimset('Display', 'none');
% 
% objective = @(x) -resample_optim(data, x(1), x(2), 10);
% 
% % optimize
% x = fmincon(objective, x0, [], [], [], [], x0-0.5, x0+0.5, [], options);

[M, N] = size(cspace);
phi = (0:M-1)'*(0:N-1); 

% return 
rpk = [];
cpk = [];
stat = [];

end

function val = resample_sinc(CC, row, col)
% function val = resample_sinc(CC, row, col)
%
% Resample the input matrix using a lanczos (windowed sinc) interpolation
% kernel.
%
% Arguments:
%
% %

end