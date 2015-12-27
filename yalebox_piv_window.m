function [win, pos, frac_data, r_centroid, c_centroid] = ...
                yalebox_piv_window(img, r_center, c_center, len)
% function [win, pos, frac_data, r_centroid, c_centroid] = ...
%                 yalebox_piv_window(img, r_center, c_center, len)            
%
% Extract and return a sample or interrogation window from the input image,
% padding as needed, the fraction of the window that contains data (~= 0), and
% the centroid of the data. 
% 
% Regions without data are filled with random white noise to reduce thier
% contribution to the normalized cross correlation.
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
%       of all the pixels)
%
%   len = Scalar, integer, length of the window
%
%   win = 2D matrix, double, subset of img, possibly with zero padding
%
%   pos = Vector, length == 4, position vector for 'win' in the format
%       [left, bottom, width, height] in pixel coordinates
%
%   frac_data = Scalar, range [0, 1], fraction of the window that contains data
%       (i.e. sand)
%
%   r_centroid, c_centroid = Scalar, centroid of the data (i.e. sand) in the
%       window
% %

% get window limits, may lie outside the image domain
hlen = (len-1)/2;
r0 = floor(r_center-hlen);
r1 =  ceil(r_center+hlen); 
c0 = floor(c_center-hlen);
c1 =  ceil(c_center+hlen);

% generate position vector for output
pos = [c0, r0, c1-c0, r1-r0];

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
         
% replace no-data pixels with random noise too
no_data = win==0;
% rand_win = rand(size(win));
% win(no_data) = rand_win(no_data);

% compute fraction of the window that has data
frac_data = 1-sum(no_data(:))/numel(win);

if nargout >= 4
    % compute centroid of the data in the window
    [r_data, c_data] = find(~no_data);
    n_data = length(r_data);
    r_centroid = sum(r_data)/n_data;
    c_centroid = sum(c_data)/n_data;    
    
%     fprintf('r_centroid = %f, c_centroid = %f\n', r_centroid, c_centroid);
%     imagesc(win); hold on; plot(c_centroid, r_centroid, 'xk'); hold off
    
    % convert centroid from local coordinates to full (parent) matrix coordinates
    r_centroid = r_centroid+pos(2)-1;
    c_centroid = c_centroid+pos(1)-1;
    
    
    
end