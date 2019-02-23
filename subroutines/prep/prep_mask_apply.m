function mask = prep_mask_apply(rgb, model, show)
% function mask = prep_mask_apply(rgb, model, show)
% 
% Apply trained classifier and morphological clean-up to create a
% sand/other mask array
% 
% Arguments:
%   rgb: 3D matrix, RGB 24-bit image
% 
%   model: ML model class, trained classifier
% 
%   show: set true to display resulting masked image for inspection, or
%       (default) false to skip the plot
% 
% Returns:
%   mask: 2D matrix, 1 where sand, 0 where other
% % 

% TODO: add log messages
% TODO: post process the mask with morph filters (like before)

% set defaults
if nargin < 3; show = false; end

% sanity checks
narginchk(2, 3);
validateattributes(rgb, {'numeric'}, {'3d'});
% TODO: validate model arg
validatestring(show, {'logical'}, {'scalar'});

% predict
X = prep_mask_features(rgb);
labels = predict(model, X);
if iscell(labels)
    % random forest models return silly cell arrays, conver to numeric
    labels = double(cell2mat(labels));
    labels(labels == double('1')) = true;
    labels(labels == double('0')) = false;
end
labels = reshape(labels, size(rgb, 1), size(rgb, 2));

% display
if show
    hf = figure
    hf.Name = 'Masking Results';
    
    hax1 = subplot(2, 1, 1)
    imagesc(rgb, 'AlphaData', labels == 1);
    axis equal tight
    hax1.Color = [1.0, 0.25, 0.25];
    hax1.YDir = 'normal';
    
    hax2 = subplot(2, 1, 2)
    imagesc(rgb, 'AlphaData', labels == 0);
    axis equal tight
    hax2.Color = [1.0, 0.25, 0.25];
    hax2.YDir = 'normal';
    
    linkaxes([hax1, hax2], 'xy');
end