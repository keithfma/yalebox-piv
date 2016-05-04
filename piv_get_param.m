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

% set parameters
samplen = [30, 30];
sampspc = [15, 15];
intrlen = [120, 60];
npass = [1, 3];
valid_max = 2;
valid_eps = 0.1;
lowess_span_pts = 9;
spline_tension = 0.95;
min_frac_data = 0.8; 
min_frac_overlap = 0.5;

% run piv
[xx, yy, uu, vv, roi] = ...
    piv(ini, fin, ini_mask, fin_mask, xw, yw, samplen, sampspc, ...
        intrlen, npass, valid_max, valid_eps, lowess_span_pts, spline_tension, ...
        min_frac_data, min_frac_overlap, true); 

% plot results
%...get data limits
[xxg, yyg] = meshgrid(xx, yy);
mask_uv = ~isnan(uu) & ~isnan(vv);
uv_xlim = [min(xxg(mask_uv)), max(xxg(mask_uv))];
uv_ylim = [min(yyg(mask_uv)), max(yyg(mask_uv))];
%... displacement magnitude and direction, [mm/step]
figure
subplot(3,1,1);
mm = sqrt(uu.^2+vv.^2);
imagesc(xx, yy, mm*1000, 'AlphaData', ~isnan(mm)); 
colorbar;
set(gca, 'YDir', 'normal', 'XLim', uv_xlim, 'YLim', uv_ylim);
title('Magnitude');
%... x-direction displacement magnitude, [mm/step]
subplot(3,1,2);
imagesc(xx, yy, uu*1000, 'AlphaData', ~isnan(uu)); 
colormap(gca, flipud(colormap)); % flow is in negative x direction
colorbar;
set(gca, 'YDir', 'normal', 'XLim', uv_xlim, 'YLim', uv_ylim);
title('U');
%... y-direction displacement magnitude, [mm/step]
subplot(3,1,3);
imagesc(xx, yy, vv*1000, 'AlphaData', ~isnan(vv)); 
colorbar;
set(gca, 'YDir', 'normal', 'XLim', uv_xlim, 'YLim', uv_ylim);
title('V');

%% Save parameters to mat file