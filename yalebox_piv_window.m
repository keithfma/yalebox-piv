function [win, pos, frac_data] = yalebox_piv_window(img, rcnt, ccnt, len)
% function [win, pos, frac_data] = yalebox_piv_window(img, rcnt, ccnt, len)
%
% Extract and return a sample or interrogation window from the input image,
% padding as needed, and the fraction of the window that contains data (~= 0).
% Regions without data are filled with random white noise to reduce thier
% contribution to the normalizae cross correlation.
%
% For sample windows: Sample grid centroids are chosen such that sample window
% edges lie at integer pixel coordinates. Rounding operations thus have no
% effect in this case.
%
% For interrogation windows: Window edges are not guaranteed to be integer pixel
% coordinates. Rounding is used to expand window extent outward to the nearest
% integer pixel coordinates. Note that this does not change the window centroid,
% since expansion is always symmetric.
%
% Arguments:
%
%   img = 2D matrix, double, range 0 to 1, normalized grayscale image from
%       the pair to be analyzed. Should be the initial image for the sample
%       window and the final image for the interrogation window.
%
%   rcnt, ccnt = Scalar, double, location of the window centroid 
%
%   len = Scalar, integer, length of the window
%
%   win = 2D matrix, double, subset of img, possibly with zero padding
%
%   pos = Vector, length == 4, position vector for 'win' in the format
%       [left, bottom, width, height] in pixel coordinates
% %

% get window limits, may lie outside the image domain
hlen = (len-1)/2;
r0 = floor(rcnt-hlen);
r1 =  ceil(rcnt+hlen); 
c0 = floor(ccnt-hlen);
c1 =  ceil(ccnt+hlen);

% generate position vector for output
pos = [c0, r0, c1-c0+1, r1-r0+1];

% get pad size, restrict window indices to valid range
pl = max(0, 1-c0);
c0 = max(1, c0);

nc = size(img, 2);
pr = max(0, c1-nc);
c1 = min(nc, c1);

pb = max(0, 1-r0);
r0 = max(1, r0);

nr = size(img, 1);
pt = max(0, r1-nr);
r1 = min(nr, r1);

% extract data and add pad
sub = img(r0:r1, c0:c1);
[snr, snc] = size(sub);
win = [zeros(pb, pl+snc+pr);
       zeros(snr, pl), sub, zeros(snr, pr);
       zeros(pt, pl+snc+pr)];  

% compute fraction of the window that has data
frac_data = sum(win(:)~=0)/numel(win);
