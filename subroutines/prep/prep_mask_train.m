function model = prep_mask_train(rgb, sand_poly, other_poly, model_type)
% function model = prep_mask_train(rgb, sand_poly, other_poly, model_type)
% 
% Train classifier to label pixels as sand (1) or other (0)
%
% Arguments:
%   rgb: 3D matrix, RGB 24-bit image
% 
%   sand_poly: 2D array, vertices of mask
%       polygons for sand-pixel training data, x-coords in row 1 and
%       y-coords in row 2, polygons separated by NaN
% 
%   other_poly: 2D array, vertices of mask
%       polygons for sand-pixel training data, x-coords in row 1 and
%       y-coords in row 2, polygons separated by NaN
% 
%   model_type: string, selects model from a few options: 'tree' uses a
%       simple decision tree (relatively fast for development), and 'forest'
%       uses a slower, better random forest model.
% 
%  Returns:
%   model: ML model class, trained classifier
% % 

% TODO: add log messages
% TODO: post process the mask with morph filters (like before)

% set defaults
if nargin < 4; model_type = 'tree'; end

% sanity check
validateattributes(rgb, {'numeric'}, {'3d'});
validateattributes(sand_poly, {'numeric'}, {'2d'});
validateattributes(other_poly, {'numeric'}, {'2d'});
validatestring(model_type, {'forest', 'tree'});

% get training set indices
% TODO: update upstream function so that masks are not inverted
sand_mask = ~prep_mask_manual(rgb, sand_poly);
sand_idx = find(sand_mask);
other_mask = ~prep_mask_manual(rgb, other_poly);
other_idx = find(other_mask);

% balance classes by random duplication of the smaller class
if length(sand_idx) < length(other_idx)
    num_extra = length(other_idx) - length(sand_idx);
    extra_idx = randsample(sand_idx, num_extra, true);
    sand_idx = [sand_idx; extra_idx];
else
    num_extra = length(sand_idx) - length(other_idx);
    extra_idx = randsample(other_idx, num_extra, true);
    other_idx = [other_idx; extra_idx];    
end

% extract training set
X = prep_mask_features(rgb);
X_train = X([sand_idx; other_idx], :);
Y_train = [ones(size(sand_idx)); zeros(size(other_idx))];

% fit model
switch model_type
    case 'forest'
        model = TreeBagger(10, X_train, Y_train, ...
                           'InBagFraction', 0.80, ...
                           'Cost', [0, 1; 1, 0], ...
                           'NumPrint', 1); %  , ...
                         % 'Options', statset('UseParallel', true)); 
    case 'tree'
        model = fitctree(X_train, Y_train, ...
                         'Cost', [0, 1; 1, 0]);
    
    otherwise
        error('Bad value for input argument "model_type"');
end