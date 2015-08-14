function [s] = nansmooth2(f,rng)
% function [s] = nansmooth2(f,rng)
%
% Moving average filter ignoring NaNs.  Simple implementation averages as
% many of the neighbors within rng as possible. Useful to smooth within
% irregular regions (bounded by NaNs).
%
% f = 2D Matrix, input matrix to be smoothed
%
% rng = Scalar.  Number of neighbors in each direction to include in
% averaging filter.  For example, n=1 will average a square of 8 neighbors
% and 1 center element (where available), n=2 will average 15 nieghbors and
% one central cell, etc.  Default is 1.
%
% Keith Ma, Yale University 2012

if nargin<2
    rng = 1;
end

[nr nc] = size(f);
s = nan(nr,nc);

for r = 1:nr
    for c = 1:nc
        
        if isnan(f(r,c)); % skip if the central value is a NaN
            s(r,c) = NaN;            
        else            
            sub = f(max(r-rng,1):min(r+rng,nr),max(c-rng,1):min(c+rng,nc)); % center and 8-neighbors (if available)
            sub = sub(~isnan(sub)); % drop NaNs
            s(r,c) = mean(sub); % compute mean
        end
        
    end
end
   