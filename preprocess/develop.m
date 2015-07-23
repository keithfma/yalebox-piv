% Script for developing and debugging preprocessing tools.

% NOTES:

%% define input variables

img_file = 'test_data/test_b_04.png';
entropy_radius = 9;
open_radius = 36; 
bw_threshold = 0.75;
loess_frac = 0.25;

img = imread(img_file);
img = double(img);
img = img./255;

%% test

mask = correct_mask(img, entropy_radius, open_radius, bw_threshold);
img(repmat(~mask, [1, 1, 3])) = 0;
