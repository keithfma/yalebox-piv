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
narginchk(3, 4)

validateattributes(hsv, {'double'}, {'3d', '>=', 0, '<=', 1}, ...
    'yalebox_prep_intensity', 'hsv');

validateattributes(mask, {'logical'}, {'2d', 'size', [size(hsv,1), size(hsv, 2)], ...
    'yalebox_prep_intensity', mask);

validateattributes(num_tiles, {'numeric'}, {'integer', 'positive', 'numel', 2}, ...
    'yalebox_prep_intensity', num_tiles);

if nargin == 3; show = false; end
validateattributes(show, {'numeric', 'logical'}, {'scalar'}, ...
    'yalebox_prep_world_coord', 'show');

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