function val = interp_point_fourier2(f, r0, c0)
%
% Interpolate a single point in a real 2D matrix using the fourier
% transform and the shift theorum.

% The size of the matrix.
[N, M] = size(f);

% FFT of our possibly padded input signal.
F = fft2(f);

% The mathsy bit. The floors take care of odd-length signals.
r_shift = exp(-1i * 2 * pi * r0 * [0:floor(N/2)-1 floor(-N/2):-1]' / N);
c_shift = exp(-1i * 2 * pi * c0 * [0:floor(M/2)-1 floor(-M/2):-1]  / M);

% Force conjugate symmetry. Otherwise this frequency component has no
% corresponding negative frequency to cancel out its imaginary part.
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

% There should be no imaginary component (for real input
% signals) but due to numerical effects some remnants remain.
val = real(val);