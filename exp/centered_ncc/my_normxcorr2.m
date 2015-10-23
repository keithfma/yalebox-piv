function [C] = my_normxcorr2(X, H)
% function [C] = my_normxcorr2(X, H)
%
% Pure MATLAB implementation of the space-domain normalized cross correlation.
% Based on the R2015b xcorr2 and normxcorr2 documentation, so I use the same
% variable names where possible.

% initialize
[M, N] = size(X);
[P, Q] = size(H);

mm = 1:M;
nn = 1:N;

C = zeros(M+P-1, N+Q-1);

H_pad = padarray(H, [M-1, N-1], 0, 'both');

X = X(:);
X_hat = X-mean(X);
X_var = var(X);

for kk = -(P-1):(M-1);
    for ll = -(Q-1):(N-1);
        H_shift = H_pad(mm-kk+M-1, nn-ll+N-1);
        H_shift = H_shift(:);
        H_shift_hat = H_shift-mean(H_shift);
        H_shift_var = var(H_shift);
        C(kk+P,ll+Q) = sum(X_hat.*H_shift_hat)/sqrt(X_var*H_shift_var+eps)/(M*N-1); 
        if isnan( C(kk+P,ll+Q)); keyboard; end
    end
end

% invert dimensions to match standard NCC
C = C(end:-1:1, end:-1:1);