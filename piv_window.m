function [win, r0, c0, r_centroid, c_centroid, frac_data] = piv_window(img, r_center, c_center, len)
% function [win, frac_data, r0, c0, r_centroid, c_centroid] = ...
%     piv_window(img, r_center, c_center, len)
%
% Extract and return a sample or interrogation window from the input image,
% padding as needed, the fraction of the window that contains data (~= 0), and
% the centroid of the data. 
% 
% Regions without data are filled with 0's to facilitate masking during
% normalized cross correlation.
%
% For interrogation windows: Window edges are not guaranteed to be integer pixel
% coordinates. Rounding is used to expand window extent outward to the nearest
% integer pixel coordinates. 
%
% For sample windows: Sample grid centroids are chosen such that sample window
% edges lie at integer pixel coordinates. Rounding operations are still
% performed, but have no effect in this case.
%
% Arguments:
%
%   img = 2D matrix, double, range 0 to 1, normalized grayscale image from
%       the pair to be analyzed. Should be the initial image for the sample
%       window and the final image for the interrogation window.
%
%   r_center, c_center = Scalar, double, location of the window center (centroid
%       of all the pixels, including those outside the ROI)
%
%   len = Scalar, integer, length of the window
%
%   win = 2D matrix, double, subset of img, possibly with zero padding
%
%   r0, c0 = Scalar, row and column indices of the bottom left corner of the
%       window, including any rounding
% 
%   r_centroid, c_centroid = Scalar, centroid of the data (i.e. sand) in the
%       window
%
%   frac_data = Scalar, range [0, 1], fraction of the window that contains data
%       (i.e. sand)
% %

% get window limits, may lie outside the image domain
hlen = (len-1)/2;
r0 = floor(r_center-hlen);
r1 =  ceil(r_center+hlen); 
c0 = floor(c_center-hlen);
c1 =  ceil(c_center+hlen);

% get pad size, restrict window indices to valid range
pl = max(0, 1-c0);
c0 = max(1, c0);

nc = size(img, 2);
pr = max(0, c1-nc);
c1 = min(nc, c1);

pb = max(0, 1-r0); % incorrect for negative r0
r0 = max(1, r0);

nr = size(img, 1);
pt = max(0, r1-nr);
r1 = min(nr, r1);

% extract data and add (no data) pad
sub = img(r0:r1, c0:c1);
[snr, snc] = size(sub);
win = [zeros(pb, pl+snc+pr);
    zeros(snr, pl), sub, zeros(snr, pr);
    zeros(pt, pl+snc+pr)];

if nargout > 3
    % compute centroid of the data in the window
    no_data = win==0;
    [r_data, c_data] = find(~no_data);
    n_data = length(r_data);
    r_centroid = sum(r_data)/n_data;
    c_centroid = sum(c_data)/n_data;    
    
    % convert centroid from local coordinates to full (parent) matrix coordinates
    r_centroid = r_centroid+r0-1;
    c_centroid = c_centroid+c0-1;
    
end

if nargout == 6
    % compute fraction of the window that has data

    frac_data = 1-sum(no_data(:))/numel(win);
end