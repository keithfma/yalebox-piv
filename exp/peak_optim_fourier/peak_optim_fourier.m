function [rpk, cpk, ok] = peak_optim_fourier(fg)
%
% Non-linear optimization to find interpolated peak in the xcor plane. See
% [1] for details. Using varaible names from [1] for now.
%
% References:
% 
% [1] Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
%   "Efficient subpixel image registration algorithms," Opt. Lett. 33,
%   156-158 (2008)
%
% %

FG = fft2(fg);

% debug {
subplot(1,4,1); 
imagesc(fg);
title('fg');

subplot(1,4,2); 
imagesc(real(FG)); 
title('real(FG)');

subplot(1,4,3); 
imagesc(imag(FG));
title('imag(FG)');

subplot(1,4,4); 
imagesc(ifft2(FG));
title('fg, recovered');

keyboard
% } debug
