function [rpk, cpk, ok] = peak_optim_fourier(f)
%
% Non-linear optimization to find interpolated peak in the xcor plane. See
% [1] for discussion of this approach. Implementation of the sampling
% function is modified from [2] to compute the shifted 2D FFT at one point
% only.
%
% References:
% 
% [1] Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
%   "Efficient subpixel image registration algorithms," Opt. Lett. 33,
%   156-158 (2008)
%
% [2] http://www.mathworks.com/matlabcentral/fileexchange/23440-2d-fourier-shift
%
% %

% find peak for initial guess
[rpk, cpk] = find(f == max(f(:)));

if numel(rpk) ~= 1 || numel(cpk) ~= 1
    ok = false;
    return
end

% initialize optimization
F = fft2(f);
objective = @(p) -sample_fft2(F, p(1), p(2));
x0 = [rpk, cpk];
opt = optimset(...
    'Algorithm', 'interior-point', ...
    'Display', 'off', ...
    'UseParallel', false);

% optimize
x = fmincon(objective, x0, [], [], [], [], x0-1, x0+1, [], opt); 

% return
rpk = x(1);
cpk = x(2);
ok = true;

end

function val = sample_fft2(F, r0, c0)
%
% Sample (interpolate) a single point from the fourier transform of a real
% matrix. Implementation is modified from [2].
%
% %

[N, M] = size(F);

% evaluate exp() functions for the inverse DFT at the desired point 
% ...the floors take care of odd-length signals
exp_val_r = exp(1i * 2 * pi * (r0-1) * [0:floor(N/2)-1 floor(-N/2):-1]' / N);
exp_val_c = exp(1i * 2 * pi * (c0-1) * [0:floor(M/2)-1 floor(-M/2):-1]  / M);

% force conjugate symmetry. Otherwise this frequency component has no
% corresponding negative frequency to cancel out its imaginary part.
if mod(N, 2) == 0
	exp_val_r(N/2+1) = real(exp_val_r(N/2+1));
	exp_val_c(M/2+1) = real(exp_val_c(M/2+1));
end

% complete the inverse DFT at the desired point 
Fexp = F .* (exp_val_r * exp_val_c);
val = sum(Fexp(:))/(M*N);

% There should be no imaginary component (for real input
% signals) but due to numerical effects some remnants remain.
val = real(val);

end