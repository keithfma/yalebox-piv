function img = yalebox_prep_intensity(rgb, mask, num_tiles, show)
% function img = yalebox_prep_intensity(rgb, mask, show)
%
% Apply local histogram-equalization color correction, ignoring non-sand
% pixels. This is done using a modified version of MATLAB's adapthisteq()
% tool, which is included below. See "help adapthisteq" for more details.
%
% Arguments:
%
%   rgb = 3D matrix, uint8, a 24-bit "Truecolor" RGB image, as read into
%       MATLAB with imread()
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
assert(isa(rgb, 'uint8') && size(rgb,3) == 3, ...
     'rgb is not a Truecolor (24-bit) RGB image');
assert(isa(mask, 'logical'), ...
    'mask is not of type logical ');
assert(size(mask,1) == size(rgb, 1) && size(mask,2) == size(rgb, 2), ...
    'mask and rgb dimensions do not match');
assert(numel(num_tiles) == 2 & ...
    mod(num_tiles(1), 1) == 0 & ...
    mod(num_tiles(2), 1) == 0, ...
    'num_tiles is not a two-element vector of integers.');
if nargin == 5; show = false; end
assert(numel(show) == 1 && (show == 0 || show == 1), ...
    'show is not a logical flag');

% pre-defined parameters
clip_limit = 1;
num_bins = 1e4;

% convert image to 2D intensity matrix
hsv = rgb2hsv(rgb);
val = hsv(:,:,3);

% local histogram equalization (CLAHE) ignoring non-sand pixels
img = adapthisteq(val, ...
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
    imagesc(val); 
    caxis([0,1]);
    subplot(2,1,2)
    title('equalized intensity')
    imagesc(img)
    caxis([0,1]);
    colormap(gray)
end