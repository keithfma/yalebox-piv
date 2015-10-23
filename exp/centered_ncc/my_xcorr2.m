function [C] = my_xcorr2(X, H)
% function [C] = my_xcorr2(X, H)
%
% Pure MATLAB implementation of the space-domain cross correlation. Based
% on the R2015b xcorr2 documentation, so I use the same variable names
% where possible.

% initialize
[M, N] = size(X);
[P, Q] = size(H);

mm = 1:M;
nn = 1:N;

C = zeros(M+P-1, N+Q-1);

H_pad = padarray(H, [M-1, N-1], 0, 'both');

for kk = -(P-1):(M-1);
    for ll = -(Q-1):(N-1);
        H_pad_sub = H_pad(mm-kk+M-1, nn-ll+N-1);
        C(kk+P,ll+Q) = sum(X(:).*H_pad_sub(:)); % sum() reduces numerical roundoff error   
        % C(kk+P,ll+Q) = X(:)'*H_pad_sub(:);        
    end
end