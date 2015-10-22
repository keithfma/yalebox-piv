function [C] = my_xcorr2(X, H)
% function [C] = my_xcorr2(X, H)
%
% Pure MATLAB implementation of the space-domain cross correlation. Based
% on the R2015b xcorr2 documentation, so I use the same variable names
% where possible.

% input matrix dimensions
[M, N] = size(X);
[P, Q] = size(H);

% pad H
Hp = padarray(H, [M-1, N-1], 0, 'both');

% compute 
C = zeros(M+P-1, N+Q-1);
for kk = -(P-1):(M-1);
    for ll = -(Q-1):(N-1);
        for mm = 1:M
            for nn = 1:N                 
                C(kk+P,ll+Q) = C(kk+P,ll+Q)+X(mm,nn)*Hp(mm-kk+M-1, nn-ll+N-1);
            end
        end
    end
end

% keyboard