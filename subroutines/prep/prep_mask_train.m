function model = prep_mask_train(rgb, labels_poly, model_type)
% function model = prep_mask_train(rgb, labels_poly, model_type)
% 
% Train classifier to label pixels as sand (1) or other (0)
%
% Arguments:
%   rgb: 3D matrix, RGB 24-bit image
% 
%   labels_poly: TODO
% 
%   model_type: string, select mask model from a few options: 'tree' uses a
%       simple decision tree (relatively fast for development), and 'forest'
%       uses a slower, better random forest model.
% 
%  Returns:
%   model: ML model class, trained classifier
% % 

% TODO: add log messages
% TODO: post process the mask with morph filters (like before)

% set defaults
if nargin < 3; model_type = 'tree'; end

% sanity check
validateattributes(rgb, {'numeric'}, {'3d'});
% TODO: check labels_poly
validatestring(model_type, {'forest', 'tree', 'naiive_bayes'});

% get training set labels as linear indices into a band of rgb
[labels, ~] = prep_mask_labels(rgb, labels_poly);

% count number of labels in each class
unique_labels = unique(labels(labels > 0)); % skip no-label ID
count_labels = zeros(length(unique_labels), 1);
for ii = 1:length(unique_labels)
    count_labels(ii) = sum(labels(:) == unique_labels(ii));
end

% balance classes by random resampling to desired size
num_classes = length(unique_labels);
num_balanced = 1e4;
label_idx = zeros(num_balanced, num_classes, 'int64');
for ii = 1:num_classes
   label_idx(:, ii) = randsample(find(labels == unique_labels(ii)), num_balanced, true);
end

% extract training set
X = prep_mask_features(rgb);
X_train = nan(num_classes*num_balanced, size(X, 2));
Y_train = nan(num_classes*num_balanced, 1);
for kk=1:num_classes
    i1 = kk*num_balanced;
    i0 = 1 + i1 - num_balanced;
    X_train(i0:i1, :) = X(label_idx(:, kk), :);
    Y_train(i0:i1) = kk;
end

% fit model
% TODO: fiddle with model parameters
switch model_type
    case 'forest'
        model = TreeBagger(20, X_train, Y_train, ...
                           'InBagFraction', 0.80, ...
                           'NumPrint', 1); %  , ...
                         % 'Options', statset('UseParallel', true)); 
    case 'tree'
        model = fitctree(X_train, Y_train);

    case 'naiive_bayes'
        model = fitcnb(X_train, Y_train);

                     
    otherwise
        error('Bad value for input argument "model_type"');
end