function img = yalebox_prep_intensity(hsv, mask, num_tiles, show)
% function img = yalebox_prep_intensity(hsv, mask, num_tiles, show)
%
% Apply local histogram-equalization color correction, ignoring non-sand
% pixels. This is done using a modified version of MATLAB's adapthisteq()
% tool, which is included below. See "help adapthisteq" for more details.
%
% Arguments:
%
%   hsv = 3D matrix, double, Color image in HSV colorspace.
%
%   mask = 2D matrix, double, TRUE where there is sand, FALSE elsewhere
%
%   num_tiles = 2-element vector, double, parameter to adapthisteq. From
%       the MATLAB help: "Two-element vector of positive integers specifying
%       the number of tiles by row and column, [M N]. Both M and N must be at
%       least 2. The total number of tiles is equal to M*N."
%   
%   show = Optional, scalar, logical flag, set to True to plot the original and
%       equalized intensity images, default = false
%
%   img = 2D matrix, size(mask), double, normalized intensity image derived from
%     rgb, histogram is approximately uniform
%
% Keith Ma, July 2015

% check for sane inputs
assert(isa(hsv, 'double') && size(hsv,3) == 3, ...
    'hsv is not a HSV image');
assert(isa(mask, 'logical'), ...
    'mask is not of type logical ');
assert(size(mask,1) == size(hsv, 1) && size(mask,2) == size(hsv, 2), ...
    'mask and rgb dimensions do not match');
assert(numel(num_tiles) == 2 & ...
    mod(num_tiles(1), 1) == 0 & ...
    mod(num_tiles(2), 1) == 0, ...
    'num_tiles is not a two-element vector of integers.');
if nargin == 3; show = false; end
assert(numel(show) == 1 && (show == 0 || show == 1), ...
    'show is not a logical flag');

% pre-defined parameters
clip_limit = 1;
num_bins = 1e4;

% local histogram equalization (CLAHE) ignoring non-sand pixels
img = masked_adapthisteq(hsv(:,:,3), mask, ...
        'NumTiles', num_tiles, ...
        'NBins', num_bins, ...
        'ClipLimit', clip_limit, ...
        'Range', 'full');
    
% add a tiny offset so that no sand pixels are exactly zero
tiny = 1e-5;
img(mask & img==0) = img(mask & img==0)+tiny;

% (optional) show original and equalized intensity image
if show
    figure()
    subplot(2,1,1)
    title('original intensity')
    imagesc(hsv(:,:,3)); 
    caxis([0,1]);
    subplot(2,1,2)
    title('equalized intensity')
    imagesc(img)
    caxis([0,1]);
    colormap(gray)
end