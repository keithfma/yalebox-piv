function eql = prep_intensity(rgb, mask, eql_len, show, verbose)
% function eql = prep_intensity(rgb, mask, eql_len, show, verbose)
%
% Convert masked color image to "equalized" grayscale image. 
%
% Equalization transforms rgb to grayscale with approximatley uniform pixel
% intensity distribution. A simple adaptive (local window) method is used, in
% which the value of each pixel is replaced with its percentile in a local
% window. This is equivalent to the normal method of using a CDF as a transform
% function, but more efficient, since it is not necessary to compute the CDF to
% convert just one pixel.
%
% Arguments:
%
%   rgb = 3D matrix, RGB color image
%
%   mask = 2D matrix, double, TRUE where there is sand, FALSE elsewhere
%
%   eql_len = Scalar, integer, odd. Side length (in pixels) for the local
%       neighborhood used to compute the transform for each pixel.
%
%   show = Optional, scalar, logical flag, set to True to plot the original and
%       equalized intensity images and some simple comparison statistics,
%       default = false
%
%   verbose = Optional, logical flag, display verbose messages (1) or don't (0)
%
%   eql = 2D matrix, size(mask), double, normalized intensity image with uniform
%       distribution in the range [0, 1].
%
% %

% set defaults
narginchk(3, 5)
if nargin < 4; show = false; end
if nargin < 5; verbose = false; end

% check for sane inputs
validateattributes(rgb, {'numeric'}, {'3d'});
validateattributes(mask, {'logical'}, {'2d', 'size', [size(rgb,1), size(rgb, 2)]});
validateattributes(eql_len, {'numeric'}, {'integer', 'positive', 'odd'});
validateattributes(show, {'numeric', 'logical'}, {'scalar'});
validateattributes(verbose, {'numeric', 'logical'}, {'scalar'});

% convert rgb to masked grayscale
img = rgb2hsv(rgb);
img = img(:,:,3); % grayscale layer is  HSV "value"
img(~mask) = 0;

% get constants, preallocate variables
[nr, nc] = size(img);
eql = zeros(nr, nc);
nhalfwin = floor(eql_len/2);

% (optional) report basic stats on intial image
if verbose
    fprintf('%s: input image: min = %.1e, mean = %.1e, max = %.1e. std = %.1e\n', ...
        mfilename, min(img(mask)), mean(img(mask)), max(img(mask)), std(img(mask)));
end

% local histogram equalization 
for i = 1:nr 
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


% (optional) report basic stats on final image
if verbose
    fprintf('%s: output image: min = %.1e, mean = %.1e, max = %.1e. std = %.1e\n', ...
        mfilename, min(eql(mask)), mean(eql(mask)), max(eql(mask)), std(eql(mask)));
end


