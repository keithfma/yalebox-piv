function eql = prep_intensity(img, mask, nwin, show)
% function eql = prep_intensity(img, mask, nwin, show)
%
% Convert masked color image to "equalized" grayscale image. 
%
% Equalization transforms pixel intensity so that the image histogram is
% approximately uniform. A simple adaptive (local window) method is used, in
% which the value of each pixel is replaced with its percentile in a local
% window. This is equivalent to the normal method of using a CDF as a transform
% function, but more efficient, since it is not necessary to compute the CDF to
% convert just one pixel.
%
% Arguments:
%
%   img = 2D matrix, double, grayscale image in the range [0, 1].
%
%   mask = 2D matrix, double, TRUE where there is sand, FALSE elsewhere
%
%   nwin = Scalar, integer, odd. Side length (in pixels) for the local
%       neighborhood used to compute the transform for each pixel.
%
%   show = Optional, scalar, logical flag, set to True to plot the original and
%       equalized intensity images and some simple comparison statistics,
%       default = false
%
%   eql = 2D matrix, size(mask), double, normalized intensity image with uniform
%       distribution in the range [0, 1].
%
% %

% set defaults
narginchk(3, 4)
if nargin == 3; show = false; end

% check for sane inputs
validateattributes(img, {'double'}, {'2d', '>=', 0, '<=', 1});
validateattributes(mask, {'logical'}, {'2d', 'size', [size(img,1), size(img, 2)]});
validateattributes(nwin, {'numeric'}, {'integer', 'positive', 'odd'});
validateattributes(show, {'numeric', 'logical'}, {'scalar'});

% get constants, preallocate variables
[nr, nc] = size(img);
eql = zeros(nr, nc);
nhalfwin = floor(nwin/2);

% apply mask
img(~mask) = 0;

% local histogram equalization 
parfor i = 1:nr 
    for j = 1:nc        
        % skip pixels outside the mask
        if mask(i,j) == 0
            continue
        end
        
        % extract window from image and mask
        i0 = max(1 , i-nhalfwin);
        i1 = min(nr, i+nhalfwin);
        j0 = max(1 , j-nhalfwin);
        j1 = min(nc, j+nhalfwin);
        img_win = img(i0:i1, j0:j1);
        mask_win = mask(i0:i1, j0:j1);
        
        % compute transform using percentile of the center point among its neighbors 
        nbr = img_win(mask_win);
        eql(i,j) = (sum(nbr<img(i,j)) + 0.5*sum(nbr==img(i,j)))/numel(nbr);        
    end
end

% add a tiny offset so that no sand pixels are exactly zero
tiny = 1e-5;
eql(mask & eql==0) = eql(mask & eql==0)+tiny;

% (optional) show original and equalized intensity image
if show
    
    % plot original grayscale and equalized image
    figure    
    subplot(2,1,1)
    imagesc(img);
    caxis([0,1]);
    axis off
    title('original grayscale');    
    subplot(2,1,2)
    imagesc(eql);
    caxis([0,1]);
    axis off
    title('equalized grayscale');

    % plot empirical CDF and histogram for original and equalized images
    figure
    nbins = 100;
    
    subplot(2,2,1)    
    [cdf, val] = ecdf(img(mask));
    plot(val, cdf, 'Marker', '.');
    title('original grayscale CDF');
    
    subplot(2,2,3)
    hist(img(mask), nbins);
    title('original grayscale histogram');
    
    subplot(2,2,2)    
    [cdf, val] = ecdf(eql(mask));
    plot(val, cdf, 'Marker', '.');
    title('equalized grayscale CDF');
    
    subplot(2,2,4)
    hist(eql(mask), nbins);
    title('equalized grayscale histogram');
    
end



