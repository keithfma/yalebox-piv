function mask = checkermask(val, sample)
% function mask = checkermask(val, sample)
% 
% Generate downsampled checkerboard mask for input value array, useful for 
% displaying overlays on image plots (e.g. quiver vectors)
%
% Arguments:
%   val: matrix, value array to generate mask for
%   sample: scalar integer, downsampling factor, must be even
% %

% assert(mod(sample, 2) == 0, 'sample argument must be even');

[nr, nc] = size(val);
mask = false(nr, nc);

for ii = 1:(2*sample):nr
    for jj = 1:(2*sample):nc
        mask(ii, jj) = true;
    end
end

for ii = (1 + sample):(2*sample):nr
    for jj = (1 + sample):(2*sample):nc
        mask(ii, jj) = true;
    end
end

