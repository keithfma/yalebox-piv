function [rpk, cpk, ok] = peak_optim_fourier(f)
%
% Non-linear optimization to find interpolated peak in the xcor plane. See
% [1] for discussion of this approach. Implementation of the sampling
% function is modified from [2] to compute the shifted 2D FFT at once point
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


function val = interp_pt_fft2(f, r0, c0)
%
% Interpolate a single point in a 2D matrix using the fourier transform and
% the shift theorum.

% The size of the matrix.
[N, M] = size(f);

% FFT of our possibly padded input signal.
F = fft2(f);

% The mathsy bit. The floors take care of odd-length signals.
x_shift = exp(-1i * 2 * pi * delta(1) * [0:floor(N/2)-1 floor(-N/2):-1]' / N);
y_shift = exp(-1i * 2 * pi * delta(2) * [0:floor(M/2)-1 floor(-M/2):-1] / M);

% Force conjugate symmetry. Otherwise this frequency component has no
% corresponding negative frequency to cancel out its imaginary part.
if mod(N, 2) == 0
	x_shift(N/2+1) = real(x_shift(N/2+1));
end 
if mod(M, 2) == 0
	y_shift(M/2+1) = real(y_shift(M/2+1));
end


Fshift = F .* (x_shift * y_shift);

% Invert the FFT.