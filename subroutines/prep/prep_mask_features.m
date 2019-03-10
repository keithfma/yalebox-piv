function features = prep_mask_features(rgb, segments, show)
% function features = prep_mask_features(rgb, segments)
% 
% Create feature matrix (rows as observations, columns as features) from
% the input experiment image
% 
% Arguments:
%   rgb: 3D matrix, RGB 24-bit image
%   segments: image segments as integer labels, same footprint as rgb,
%       output from prep_segments
%   show: set true to display plots of each computed feature for smell test
% 
% Returns:
%   features: 2D matrix, rows as observations (pixels)
% % 

% set defaults
narginchk(2, 3);
if nargin < 3; show = false; end

% sanity check
validateattributes(rgb, {'numeric'}, {'3d'});
validateattributes(segments, {'numeric'}, {'2d'});
unique_segments = unique(segments(:));
assert(all(diff(unique_segments) == 1), 'segments must be a sequential list of integers');

% generate additional input layers
fprintf('%s: generate input layers\n', mfilename);
hsv = rgb2hsv(rgb);
kernel = strel('disk', 10).Neighborhood;
v_entropy = entropyfilt(hsv(:, :, 3), kernel);

% define input layers
layer_labels = {...
    'red', ...
    'green', ...
    'blue', ...
    'hue', ...
    'saturation', ...
    'value', ...
    'value entropy'};
num_pix = size(rgb, 1)*size(rgb, 2);
layers = {...
    reshape(rgb(:,:,1), num_pix, 1), ...
    reshape(rgb(:,:,2), num_pix, 1), ...
    reshape(rgb(:,:,3), num_pix, 1), ...
    reshape(hsv(:,:,1), num_pix, 1), ...
    reshape(hsv(:,:,2), num_pix, 1), ...
    reshape(hsv(:,:,3), num_pix, 1), ...
    reshape(v_entropy, num_pix, 1)};

% more sanity checks
assert(length(layer_labels) == length(layers), "Labels and layers don't match");
num_layers = length(layers);

% compute features and gather up labels
fprintf('%s: compute features per image segment\n', mfilename);
feature_labels = cell(0);
features = cell(1, num_layers);
for ii = 1:num_layers
    [this_labels, this_features] = compute_features(...
        layer_labels{ii}, double(layers{ii}), segments(:));
    num_labels = length(this_labels);
    feature_labels(end+1:end+num_labels) = this_labels(:);
    features{ii} = this_features;
end
features = cell2mat(features);

% optional display
if show
    fprintf('%s: display features, pausing between each\n', mfilename);
    figure('Units', 'Normalized', 'Outerposition', [0 0 1 1]);
    for ii = 1:length(feature_labels)
        this_features = features(:, ii);
        imagesc(this_features(segments));
        set(gca, 'YDir', 'Normal');
        colorbar
        axis equal tight
        title(feature_labels{ii});
        pause                
    end
end


function [labels, features] = compute_features(layer_name, layer_data, segment_data)
%
% Compute features, returning a row vector of length num_features
% %

func_names = {'pct5', 'pct25' 'pct50', 'pct75', 'pct95'};
func_handle = @(x) [mean(x), std(x), prctile(x, [5, 25, 50, 75, 95])];

labels = cell(1, length(func_names));
for ii = 1:length(labels)
    labels{ii} = sprintf('%s - %s', layer_name, func_names{ii});
end

features = splitapply(func_handle, layer_data, segment_data);
