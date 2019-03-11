function model = prep_mask_train(features, labels, model_type)
% function model = prep_mask_train(features, labels, model_type)
% 
% Train classifier to label pixels as sand (1) or other (2)
%
% Arguments:
%   features: TODO
% 
%   labels: TODO
% 
%   model_type: string, select mask model from a few options: 'tree' uses a
%       simple decision tree (relatively fast for development), and 'forest'
%       uses a slower, better random forest model.
% 
%  Returns:
%   model: ML model class, trained classifier
% % 

% TODO: add log messages

% set defaults
narginchk(2, 3);
if nargin < 3; model_type = 'tree'; end

% sanity check
% TODO: complete sanity checks\
% TODO: check that labels contains only 0, 1, 2
validatestring(model_type, {'forest', 'tree', 'naiive_bayes'});

% drop no-class labels
has_class_idx = find(labels ~= 0);
features = features(has_class_idx, :);
labels = labels(has_class_idx);

% count number of labels in each class
count_labels = [0, 0];
for ii = 1:2
    count_labels(ii) = sum(labels(:) == ii);
end

% balance classes by random resampling to match larger class
if count_labels(1) < count_labels(2)
    num_resample = count_labels(2) - count_labels(1);
    fprintf('%s: balance classes by adding %i resampled sand features\n', mfilename, num_resample);
    sample_idx = randsample(find(labels == 1), num_resample, true);
elseif count_labels(2) < count_labels(1)
    num_resample = count_labels(1) - count_labels(2);
    fprintf('%s: balance classes by adding %i resampled non-sand features\n', mfilename, num_resample);
    sample_idx = randsample(find(labels == 2), num_resample, true);
else
    sample_idx = [];
end
labels = [labels; labels(sample_idx)];
features = [features; features(sample_idx, :)];    

% fit model
% TODO: fiddle with model parameters
switch model_type
    case 'forest'
        fprintf('%s: train random forest model\n', mfilename);
        model = TreeBagger(20, features, labels, ...
                           'InBagFraction', 0.80, ...
                           'NumPrint', 1); %  , ...
                         % 'Options', statset('UseParallel', true)); 
    case 'tree'
        fprintf('%s: train decision tree model\n', mfilename);
        model = fitctree(features, labels);

    case 'naiive_bayes'
        fprintf('%s: train naiive bayes model\n', mfilename);
        model = fitcnb(features, labels);

    otherwise
        error('Bad value for input argument "model_type"');
end