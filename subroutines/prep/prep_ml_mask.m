function mask = prep_ml_mask(rgb, sand_poly, other_poly)

% TODO: add log messages
% TODO: handle both training and prediction modes
% TODO: post process the mask with morph filters (like before)
% TODO: explore different values for the entropy filter size -- or just use
%   several kernels as separate features
% TODO: explore better models, rather than a shitty lil' decision tree
% TODO: try unequal loss, I care more about NOT including bad pixels than
%   losing good ones, morphology handles the fill for me
% TODO: try downsampling the input points so I can use more features?

% get training set indices
% note: masks are invereted for this application
sand_mask = ~prep_mask_manual(rgb, sand_poly);
other_mask = ~prep_mask_manual(rgb, other_poly);

sand_idx = find(sand_mask);
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
    
% compute feature values
hsv = rgb2hsv(rgb);

entropy_kernel = strel('disk', 5).Neighborhood;
h_e = entropyfilt(hsv(:, :, 1), entropy_kernel);
s_e = entropyfilt(hsv(:, :, 2), entropy_kernel);
v_e = entropyfilt(hsv(:, :, 3), entropy_kernel);

% build feature matrix
as_feature = @(x) double(reshape(x, [], 1));
features = [...
    as_feature(hsv(:, :, 1)), ...
    as_feature(hsv(:, :, 2)), ...
    as_feature(hsv(:, :, 3)), ...
    as_feature(h_e), ...
    as_feature(s_e), ...
    as_feature(v_e), ...
    ];

% train
X_train = features([sand_idx; other_idx], :);
Y_train = [ones(size(sand_idx)); zeros(size(other_idx))];

% % random forest, version 2
% model = TreeBagger(10, X_train, Y_train, ...
%     'Cost', [0, 1; 1, 0], ...
%     'MinLeafSize', 25); %  , ...
%     % 'Options', statset('UseParallel', true)); 

% decision tree with custom cost
% note:  pretty good now! fancy predictor selection does not help, needs
%   moph cleanup
model = fitctree(X_train, Y_train, ...
    'Cost', [0, 1; 1, 0], ...
    'MinLeafSize', 25); % , ...
    % 'PredictorSelection', 'interaction-curvature'); 

% % note: random forest (default) with 10 trees no better than decision tree
% model = TreeBagger(10, X_train, Y_train); 

% note: decent, but not there yet
% model = fitctree(X_train, Y_train);

% predict
labels = predict(model, features);
labels = reshape(labels, size(rgb, 1), size(rgb, 2));

% display
figure

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

disp('debug');