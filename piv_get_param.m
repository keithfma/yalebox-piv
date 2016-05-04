% Script. Template script for exploring and selecting image piv parameters. The
% intended wusage is to make a copy of the script for a given experiment, run it
% cell by cell, modifying the default parameters to suit the experiment
% particulars.

%% Init

image_in_file = '../yalebox-exp-erosion/data/K23_side.image.nc';
param_out_file = 'test/prep_get_parameters.mat';
ini_step = 10;

%% Read in image pair, masks, coordinate vectors

% find index of ini and fin
step = ncread(image_in_file, 'step');
ini_index = find(step == ini_step);

% read images and masks
mask_manual = ncread(image_in_file, 'mask_manual');
mask_manual = logical(mask_manual);

ini = ncread(image_in_file, 'image', [1, 1, ini_index],   [inf, inf, 1]); 
ini = double(ini);

fin = ncread(image_in_file, 'image', [1, 1, ini_index+1], [inf, inf, 1]); 
fin = double(fin);

ini_mask = ncread(image_in_file, 'mask_auto', [1, 1, ini_index],   [inf, inf, 1]); 
ini_mask = logical(ini_mask) & mask_manual;

fin_mask = ncread(image_in_file, 'mask_auto', [1, 1, ini_index+1],   [inf, inf, 1]); 
fin_mask = logical(ini_mask) & mask_manual;

xw = ncread(image_in_file, 'x');
xw = double(xw);

yw = ncread(image_in_file, 'y');
yw = double(yw);

% display images
figure 
ax = subplot(2,1,1);
imagesc(xw, yw, ini);
colormap('gray');
ax.YDir = 'normal';
title('ini');

ax = subplot(2,1,2);
imagesc(xw, yw, fin);
colormap('gray');
ax.YDir = 'normal';
title('fin');

figure 
ax = subplot(2,1,1);
imagesc(xw, yw, ini_mask);
colormap('gray');
ax.YDir = 'normal';
title('ini\_mask');

ax = subplot(2,1,2);
imagesc(xw, yw, fin_mask);
colormap('gray');
ax.YDir = 'normal';
title('fin\_mask');

%% Set parameters and run PIV

piv(ini, fin, ini_mask, fin_mask, xw, yw, samplen, sampspc, ...
        intrlen, npass, valid_max, valid_eps, lowess_span_pts, spline_tension, ...
        min_frac_data, min_frac_overlap, verbose) 


