function mask = prep_mask_auto(hsv, hue_lim, val_lim, entr_lim, entr_win, ...
                    morph_open_rad, morph_erode_rad, show, verbose)
% function mask = prep_mask_auto(hsv, hue_lim, val_lim, entr_lim, entr_win, ...
%                     morph_open_rad, morph_erode_rad, show, verbose)
%
% Create a logical mask for a color image that is TRUE where there is sand and
% FALSE elsewhere. This can be used to remove (set to 0) the background in a
% image prior to PIV analysis or other applications. Sand is identified by
% remapping to HSV colorspace, filtering and thresholding the "hue" and "value"
% bands.
%
% Arguments:
%
%   hsv = 3D matrix, double, Color image in HSV colorspace.
%
%   hue_lim = 2-element vector, double, range [0, 1]. [minimum, maximum] HSV
%     "hue" included as sand in the mask.
%
%   val_lim = 2-element vector, double, range [0,1]. [minimum, maximum] HSV
%     "value" included as sand in the mask.
%
%   entr_lim = 2-element vector, double, range [0, 1]. [minimum, maximum]
%     entropy included as sand in the mask. 
%
%   entr_win = scalar, integer, window size in pixels for entropy filter.
%
%   morph_open_rad = scalar, double, radius of disk structuring element used in
%     mophological opening filter. 
%
%   morph_erode_rad = scalar, double, radius of disk structuring element used in
%     mophological erosion filter.
%
%   show = Scalar, logical, set to 1 (true) to plot the mask bands, used to
%       facilitate the parameter selection process, default = false.
% 
%   verbose = Logical flag, display verbose messages (1) or don't (0)
%
%   mask = 2D matrix, logical, true where there is sand and false
%       elsewhere.
%
% Keith Ma

% set default values
if nargin < 7; show = false; end
if nargin < 8; verbose = false; end

% check for sane arguments, set default values
narginchk(7,9);
validateattributes(hsv, {'double'}, {'3d', '>=', 0, '<=', 1});
validateattributes(hue_lim, {'double'}, {'vector', 'numel', 2, '>=', 0, '<=', 1});
validateattributes(val_lim, {'double'}, {'vector', 'numel', 2, '>=', 0, '<=', 1});
validateattributes(entr_lim, {'double'}, {'vector', 'numel', 2, '>=', 0, '<=', 1});
validateattributes(entr_win, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(morph_open_rad, {'numeric'}, {'scalar', 'positive'});
validateattributes(morph_erode_rad, {'numeric'}, {'scalar', 'positive'});
validateattributes(show, {'numeric', 'logical'}, {'scalar'});
validateattributes(verbose, {'numeric', 'logical'}, {'scalar'});

% get hue, value, and entropy, normalized to the range [0, 1]
hue = hsv(:,:,1);
value = hsv(:,:,3);
entropy = entropyfilt(value, true(entr_win));

hue = hue-min(hue(:)); hue = hue./max(hue(:));
value = value-min(value(:)); value = value./max(value(:));
entropy = entropy-min(entropy(:)); entropy = entropy./max(entropy(:));

% threshold bands
hue_mask = hue >= hue_lim(1) & hue <= hue_lim(2);
value_mask = value >= val_lim(1) & value <= val_lim(2);
entropy_mask = entropy >= entr_lim(1) & entropy <= entr_lim(2);

% create mask
mask = hue_mask & value_mask & entropy_mask;

% fill holes - top and bottom
wall = true(1, size(mask,2));
mask = [wall; mask; wall];
mask = imfill(mask, 'holes');
mask = mask(2:end-1, :);

% fill holes - left and right
wall = true(size(mask, 1), 1);
mask = [wall, mask, wall];
mask = imfill(mask, 'holes');
mask = mask(:, 2:end-1);

% clean up edges with morphological filters
% ...remove noise (open), then trim boundary (erode)
mask = imopen(mask, strel('disk', morph_open_rad));
mask = imerode(mask, strel('disk', morph_erode_rad));

% (optional) plot to facilitate parameter selection
if show
    figure()
    colormap(gray);
    subplot(3,2,1); imagesc(hue); title('hue'); set(gca,'XTick', [], 'YTick',[])
    subplot(3,2,2); imagesc(hue_mask); title('hue mask'); set(gca,'XTick', [], 'YTick',[])
    subplot(3,2,3); imagesc(value); title('value'); set(gca,'XTick', [], 'YTick',[])
    subplot(3,2,4); imagesc(value_mask); title('value mask'); set(gca,'XTick', [], 'YTick',[])
    subplot(3,2,5); imagesc(entropy); title('entropy'); set(gca,'XTick', [], 'YTick',[])
    subplot(3,2,6); imagesc(entropy_mask); title('entropy mask'); set(gca,'XTick', [], 'YTick',[])
    
    figure()    
    colormap(gray);
    subplot(2,1,1); imagesc(hsv(:,:,3)); title('original'); set(gca,'XTick', [], 'YTick',[])
    subplot(2,1,2); imagesc(mask); title('mask'); set(gca,'XTick', [], 'YTick',[])    
end

% (optional) report percentage masked
if verbose
    pct_sand = 100*sum(mask(:))/numel(mask);
    fprintf('%s: %.1f%% sand, %.1f%% background\n', mfilename, pct_sand, 100-pct_sand);
end