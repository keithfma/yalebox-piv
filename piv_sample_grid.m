function [rr, cc, xx, yy] = piv_sample_grid(len, spc, xw, yw)
% function [rr, cc, xx, yy] = piv_sample_grid(len, spc, xw, yw)
%
% Create sample grid with the following characteristics:
% - sample window edges fall at integer positions
% - footprint of the sample windows is >= footprint of the image (all pixels
%   belong to one or more sample windows and no sample windows lie entirely
%   outside the image)
% - grid is centered with respect to the image (e.g. overhanging sample window
%   footprint is the same width on opposing sides, to within a pixel)
%
% Note that the first sample window requires 'len' rows and cols, and each
% additional window requires 'spc' rows and cols.
%
% Arguments:
%
%   len = Scalar, integer, length of sample window in pixels
%
%   spc = Scalar, integer, grid spacing in pixels
%
%   xw, yw = Vectors, double, world coordinates for the input images. Vector
%       sizes match the input images, such that size(image, 1) == length(yw) and
%       size(image, 2) == length(xw)
%
%   rr, cc = 2D matrix, integer, coordinate matrices, as constructed by meshgrid, 
%       for the sample grid in the y (a.k.a. row) and x (a.k.a. column) directions, 
%       in pixels
%
%   xx, yy = Vectors, double, world coordinate vectors corresponding to the
%       locations in rr, cc
% %

% get sample grid in pixel coordinates
sz = [length(yw), length(xw)];
footprint = len + ceil((sz - len)/spc)*spc; 
rem = footprint - sz;
hlen = (len - 1)/2; 

start = 1 - floor(rem/2) + hlen; 
stop = start + footprint - len;

rvec = start(1):spc:stop(1);
cvec = start(2):spc:stop(2);

[cc, rr] = meshgrid(cvec, rvec);

% interpolate pixel coordinates to world coordinates
xx = interp1(1:sz(2), xw, cvec, 'linear', 'extrap');
yy = interp1(1:sz(1), yw, rvec, 'linear', 'extrap');

% % debug: check grid edges {
% disp(rem(1))
% disp(1-(rvec(1)-hlen))
% disp((rvec(end)+hlen)-sz(1))
% disp(rem(2))
% disp(1-(cvec(1)-hlen))
% disp((cvec(end)+hlen)-sz(2))
% % } debug
