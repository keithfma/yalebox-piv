function model = prep_mask_train(features, labels)
% function model = prep_mask_train(features, labels)
% 
% Train random forest classifier to label pixels as sand (1) or other (2)
%
% Arguments:
%   features: TODO
% 
%   labels: TODO
% 
%  Returns:
%   model: ML model class, trained random forest classifier
% % 

% set defaults
narginchk(2, 2);

% sanity check
% TODO: complete sanity checks
% TODO: check that labels contains only 0, 1, 2

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
% note: minimal experimentation with model parameters
fprintf('%s: train random forest model\n', mfilename);
model = TreeBagger(100, features, labels, 'InBagFraction', 0.80); 
