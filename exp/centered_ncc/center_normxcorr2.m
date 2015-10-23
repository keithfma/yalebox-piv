function [xcr] = center_normxcorr2(samp, intr)
% function xcr = center_normxcorr2(samp, intr)
%
% Center the sample window once, and each shifted subimage of the interrogation
% window.

% initialize
X = samp;
H = intr;
[M, N] = size(X);
[P, Q] = size(H);

mm = 1:M;
nn = 1:N;

C = zeros(M+P-1, N+Q-1);

H_pad = padarray(H, [M-1, N-1], 0, 'both');

X = X(:);
X_mask = X~=0;
X_hat = X;
X_hat(X_mask) = X_hat(X_mask)-mean(X_hat(X_mask));
X_var = var(X_hat);

for kk = -(P-1):(M-1);
    for ll = -(Q-1):(N-1);
        H_shift = H_pad(mm-kk+M-1, nn-ll+N-1);
        H_shift = H_shift(:);
        H_shift_mask = H_shift~=0;
        H_shift_hat = H_shift;
        H_shift_hat(H_shift_mask) = H_shift(H_shift_mask)-mean(H_shift(H_shift_mask));
        H_shift_var = var(H_shift_hat);
        C(kk+P,ll+Q) = sum(X_hat.*H_shift_hat)/sqrt(X_var*H_shift_var+eps)/(M*N-1); 
        if isnan( C(kk+P,ll+Q)); keyboard; end
    end
end

% invert dimensions to match standard NCC
xcr = C(end:-1:1, end:-1:1);
