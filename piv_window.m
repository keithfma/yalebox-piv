function [win, r_win, c_win] = piv_window(img, r_cntr, c_cntr, len)
% function [win, r_win, c_win] = piv_window(img, r_cntr, c_cntr, len)
% 
% Extract and return a sample or interrogation window from the input image,
% padding as needed, and its coordinate information Regions without data are
% filled with 0's to facilitate masking during normalized cross correlation.
% Rounding is used to expand window extent outward to the nearest integer pixel
% coordinates.
%
% Arguments:
%
%   img = 2D matrix, double, range 0 to 1, normalized grayscale image from
%       the pair to be analyzed. Should be the initial image for the sample
%       window and the final image for the interrogation window.
%
%   r_cntr, c_cntr = Scalar, double, location of the window center (centroid
%       of all the pixels, including those outside the ROI)
%
%   len = Scalar, integer, desired length of the (square) window
%
%   win = 2D matrix, double, subset of img, possibly with zero padding
%
%   r_win, c_win = Vector, row- and column-coordinates for the returned window
% %

% get window limits, may lie outside the image domain
hlen = (len-1)/2;
r0 = floor(r_cntr - hlen);
r1 =  ceil(r_cntr + hlen); 
c0 = floor(c_cntr - hlen);
c1 =  ceil(c_cntr + hlen);

% create coordinate vectors (pixel / intrinsic units)
r_win = r0:r1;
c_win = c0:c1;

% get pad size, restrict window indices to valid range
pl = max(0, 1-c0);
c0 = max(1, c0);

nc = size(img, 2);
pr = max(0, c1-nc);
c1 = min(nc, c1);

pb = max(0, 1-r0); % assumes positive r0
r0 = max(1, r0);

nr = size(img, 1);
pt = max(0, r1-nr);
r1 = min(nr, r1);

% extract data and add (no data) pad
sub = img(r0:r1, c0:c1);
[snr, snc] = size(sub);
win = [zeros(pb, pl + snc + pr);
       zeros(snr, pl), sub, zeros(snr, pr);
       zeros(pt, pl + snc + pr)];