% Script. Template script for exploring and selecting image piv parameters. The
% intended wusage is to make a copy of the script for a given experiment, run it
% cell by cell, modifying the default parameters to suit the experiment
% particulars.

%% Init

image_in_file = '../yalebox-exp-erosion/data/K23_side.image.nc';
param_out_file = 'test/prep_get_parameters.mat';
ini_step = 10;

%% Read in image pair

% find index of ini and fin
step = ncread(image_in_file, 'step');
ini_index = find(step == ini_step);

% read images
ini = ncread(image_in_file, 'image', [1, 1, ini_index],   [inf, inf, 1]); 
fin = ncread(image_in_file, 'image', [1, 1, ini_index+1], [inf, inf, 1]); 

% convert to double
ini = double(ini);
fin = double(fin);

% display images
figure 
ax = subplot(2,1,1);
imagesc(ini);
colormap('gray');
ax.YDir = 'normal';
title('ini');

ax = subplot(2,1,2);
imagesc(fin);
colormap('gray');
ax.YDir = 'normal';
title('fin');



