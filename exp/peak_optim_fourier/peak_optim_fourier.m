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

ok = true;

% initial guess
[rpk, cpk] = find(f == max(f(:)));

if numel(rpk) ~= 1 || numel(cpk) ~= 1
    ok = false;
    return
end

% initialize optimization
F = fft2(f);

% optimize

end

function val = interp_point(F, r0, c0)
%
% Interpolate a single point from the fourier transform of a real matrix.
% Uses the fourier shift theorum to shift the point of interest to the
% origin, then computes the inverse fourier transform for that point only.
% Implementation is modified from [2].
%
% Here, this function is used to create an objective function for
% non-linear optimization.

[N, M] = size(f);

% compose shift matrix from vectors
r_shift = exp(-1i * 2 * pi * r0 * [0:floor(N/2)-1 floor(-N/2):-1]' / N);
c_shift = exp(-1i * 2 * pi * c0 * [0:floor(M/2)-1 floor(-M/2):-1]  / M);

% force conjugate symmetry, otherwise imaginary parts of frequency and
% corresponding negative frequency won't cancel
if mod(N, 2) == 0
	r_shift(N/2+1) = real(r_shift(N/2+1));
end 
if mod(M, 2) == 0
	c_shift(M/2+1) = real(c_shift(M/2+1));
end

% Apply shift
F = F .* (r_shift * c_shift);

% Invert the FFT at the shifted origin only
val = sum(F(:))/M/N;

% imaginary component isa numerical artefact for real inputs
val = real(val);

end