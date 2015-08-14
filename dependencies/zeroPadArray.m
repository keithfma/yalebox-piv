function [pA] = zeroPadArray(A,nt,nb,nl,nr)
% function [pA] = zeroPadArray(A,nt,nb,nl,nr)
%
% Pads a 2D array with the specified number of rows/cols of zeros.  If
% input array A has a 3rd dimension, the padding is repeated for all
% dimensions of A.  Dimensionality greater than 3 is not supported.
%
% Arguments-----------------
%
% A = 2D matrix.  Input array to be padded. 
%
% nt = Integer>=0.  Number of rows of zeros to append at the "top" of the
% array, where "top" refers to row positions < 1
%
% nb = Integer>=0.  Number of rows of zeros to append at the "bottom" of the
% array, where "bottom" refers to row positions > size(A,1)
%
% nl  = Integer >=0.  Number of columns of zeros to append at the "left" side of the
% array, where "left" refers to column positions < 1
%
% nr  = Integer >=0.  Number of columns of zeros to append at the "right" side of the
% array, where "right" refers to column positions > size(A,2)
%
% Output -------------------
%
% pA = Matrix. The desired padded array.
%
% Keith Ma, September 2012

% bail out for errors
if ndims(A)>3
    fprintf(2,'Error in zeroPadArray: ndims(A) must be <= 3\n');
    return
end

% allocate
[m,n,o] = size(A);
pA = nan(m+nt+nb,n+nl+nr,o); 

% pad
for i = 1:o
    
    pA(:,:,i) = [zeros(nt,n+nl+nr);...
        zeros(m,nl), A(:,:,i), zeros(m,nr);...
        zeros(nb,n+nl+nr)];
end