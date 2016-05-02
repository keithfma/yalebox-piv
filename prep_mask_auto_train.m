function cluster_center = prep_mask_auto_train(rgb, entropy_len, num_cluster, show, verbose)
% function mask = prep_mask_auto(rgb, entropy_len, num_cluster, cluster_center, show, verbose)
%
% Train a kmeans segmentation model to separate sand pixels from background.
% Input data layers are "hue" and "value" from the HSV colorspace, and a local
% entropy filter derived from the "value" layer. 
%
% Arguments:
%
%   rgb = 3D matrix, double, Color image in RGB colorspace.
%
%   entropy_len = scalar, integer, window size in pixels for entropy filter.
%
%   num_cluster = scalar, integer, number of cluster in kmeans cluster analysis.
%
%   show = Scalar, logical, set to 1 (true) to plot the cluster labels, used to
%       facilitate the parameter selection process, default = false.
% 
%   verbose = Logical flag, display verbose messages (1) or don't (0)
%
%   cluster_center = 2D matrix, cluster centers as produced by kmeans() on a
%   training image. Each row contains one cluster center in [hue, value,
%   entropy] coordinates. 
%
% Keith Ma

% set default values
if nargin < 4; show = false; end
if nargin < 5; verbose = false; end

% check for sane arguments, set default values
narginchk(3,5);
validateattributes(rgb, {'numeric'}, {'3d'});
validateattributes(entropy_len, {'numeric'}, {'scalar', 'integer', 'positive', 'odd'});
validateattributes(num_cluster, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(show, {'numeric', 'logical'}, {'scalar'});
validateattributes(verbose, {'numeric', 'logical'}, {'scalar'});

% get input data layers: hue, value, and entropy
if verbose
    fprintf('prep_mask_auto_train: computing data layers from rgb\n');
end

hsv = rgb2hsv(rgb);
hue = hsv(:,:,1);
value = hsv(:,:,3);
entropy = entropyfilt(value, true(entropy_len));

% fit the k-means cluster model
if verbose
    fprintf('prep_mask_auto_train: fitting kmeans cluster model\n');
end

cluster_data = [hue(:), value(:), entropy(:)];
[cluster_label, cluster_center] = kmeans(cluster_data, num_cluster, ...
    'Display', 'iter', ...
    'Distance', 'sqeuclidean', ...
    'Replicates', 3);
cluster_label = reshape(cluster_label, size(hue));

% interactively select the correct cluster label for the sand
if verbose
    fprintf('prep_mask_auto_train: interactive label selection\n');
end

figure;
imagesc(cluster_label);
colorbar;
while 1
    sand_label_cell = inputdlg('Select sand cluster label');
    try
        sand_label = str2double(sand_label_cell{1});
    end
    if round(sand_label) == sand_label && sand_label >= 1 && sand_label <= num_cluster
        break
    end
end
close(gcf);

% shift sand cluster to label == 1
cluster_center = circshift(cluster_center, -(sand_label-1), 1); 

% (optional) plot to facilitate parameter selection
if show
    figure()
    subplot(4,1,1); imagesc(hue); title('hue'); set(gca,'XTick', [], 'YTick',[])
    subplot(4,1,2); imagesc(value); title('value'); set(gca,'XTick', [], 'YTick',[])
    subplot(4,1,3); imagesc(entropy); title('entropy'); set(gca,'XTick', [], 'YTick',[])
    subplot(4,1,4); imagesc(cluster_label); title('cluster'); set(gca,'XTick', [], 'YTick',[])    
end