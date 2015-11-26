function [eql] = local_he_brute_force(im, ignore, nwin)
% function [eql] = local_he_brute_force(im, ignore, nwin)
%
% Local (adaptive) histogram equalization using brute-force approach.
% Equalization tranform is computed for each pixel based on a neighboring
% pixels.
%
% Arguments:
%
% im = 2D matrix, double. Continuous-valued image matrix.
%
% ignore = Scalar, double. Values to ignore in equalization routine (masked
%   pixels).
% 
% nwin = Scalar, integer, odd. Side length (in pixels) for the local
%   neighborhood used to compute the transform for each pixel.
%
% eql = 2D matrix, double. Equalized image matrix, with a uniform distribution
%   in the range [0,1]
%
% %

% check inputs
validateattributes(im, {'double'}, {'2d', 'real'}, mfilename, 'im');
validateattributes(ignore, {'double'}, {'scalar', 'real'}, mfilename, 'ignore');
validateattributes(nwin, {'numeric'}, {'integer', 'positive', 'odd'}, mfilename, 'nwin');

% initialize
[nr, nc] = size(im);
eql = zeros(nr, nc);
nhalfwin = floor(nwin/2);

% loop over all pixels
for i = 1:nr
    for j = 1:nc
        
        % skip pixels outside the region-of-interest (ROI)
        if im(i,j) == ignore
            continue
        end
        
        % extract pixels in window that are within the ROI
        i0 = max(1 , i-nhalfwin);
        i1 = min(nr, i+nhalfwin);
        j0 = max(1 , j-nhalfwin);
        j1 = min(nc, j+nhalfwin);
        win = im(i0:i1, j0:j1);
        win = win(win~=ignore);
        
        % compute transform using percentile        
        eql(i,j) = ( sum(win < im(i,j)) +0.5*sum(win == im(i,j)) )/numel(win);
        
        % % debug: monitor progress {
        % fprintf('%i, %i\n', i, j);
        % % } debug
        
    end
end