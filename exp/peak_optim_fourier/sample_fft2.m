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