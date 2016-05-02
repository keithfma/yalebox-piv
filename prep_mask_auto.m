function mask = prep_mask_auto(rgb, entropy_len, cluster_center, show, verbose)
% function mask = prep_mask_auto(rgb, entropy_len, cluster_center, show, verbose)
%
% Create a logical mask for a color image that is TRUE where there is sand and
% FALSE elsewhere. This can be used to remove (set to 0) the background in a
% image prior to PIV analysis or other applications. 
% 
% Sand is identified by segmentation using a previously trained kmeans cluster
% model. Results are cleaned up with some simple BW image operations.
%
% Arguments:
%
%   rgb = 3D matrix, double, Color image in RGB colorspace.
%
%   entropy_len = scalar, integer, window size in pixels for entropy filter.
%
%   cluster_center = 2D matrix, cluster centers as produced by kmeans() on a
%   training image. Each row contains one cluster center in [hue, value,
%   entropy] coordinates.
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
if nargin < 4; show = false; end
if nargin < 5; verbose = false; end

% check for sane arguments, set default values
narginchk(3,5);
validateattributes(rgb, {'numeric'}, {'3d'});
validateattributes(entropy_len, {'numeric'}, {'scalar', 'integer', 'positive', 'odd'});
validateattributes(cluster_center, {'numeric'}, {'2d', 'ncols', 3});
validateattributes(show, {'numeric', 'logical'}, {'scalar'});
validateattributes(verbose, {'numeric', 'logical'}, {'scalar'});

% get input data layers: hue, value, and entropy
if verbose
    fprintf('prep_mask_auto: computing data layers from rgb\n');
end

hsv = rgb2hsv(rgb);
hue = hsv(:,:,1);
value = hsv(:,:,3);
entropy = entropyfilt(value, true(entropy_len));

% apply kmeans cluster model to segment image 
if verbose
    fprintf('prep_mask_auto: segmenting image\n');
end

cluster_data = [hue(:), value(:), entropy(:)];
warning('off', 'stats:kmeans:FailedToConverge');
cluster_label = kmeans(cluster_data, size(cluster_center, 1), ...
    'MaxIter', 1,...
    'Start', cluster_center);
warning( 'on', 'stats:kmeans:FailedToConverge');
cluster_label = reshape(cluster_label, size(hue));

% create mask
mask = cluster_label == 1;

% cleanup mask
if verbose
    fprintf('prep_mask_auto_train: cleaning up mask\n');
end

% fill holes, wall off left, right and bottom
wall_lr = true(size(mask, 1), 1);
mask = [wall_lr, mask, wall_lr];
wall_b = true(1, size(mask,2));
mask = [wall_b; mask];
mask = imfill(mask, 'holes');
mask = mask(2:end, 2:end-1);

% extract largest connected object
object_label = bwlabel(mask);
largest_object = mode(object_label(object_label>0));
mask = object_label == largest_object;

% remove upper "halo"
disk = strel('disk', ceil(entropy_len/2)+1);
mask = imerode(mask, disk);

% (optional) plot to check results
if show
    figure()
    subplot(3,1,1); imagesc(value); title('value'); set(gca,'XTick', [], 'YTick',[])    
    subplot(3,1,2); imagesc(cluster_label); title('cluster'); set(gca,'XTick', [], 'YTick',[])    
    subplot(3,1,3); imagesc(mask); title('mask'); set(gca,'XTick', [], 'YTick',[])    
end