% Script. Explore and save piv parameters for experiment. Template to be copied
% and modified for specific experiments.

%% Init

image_in_file = '../data/IMAGE_FILE'; %.nc
param_out_file = 'PARAM_FILE'; %.mat
ini_step = 350;
gap = 1;

%% Read in image pair, masks, coordinate vectors

% find index of ini and fin
step = ncread(image_in_file, 'step');
ini_index = find(step == ini_step);

% read images and masks
mask_manual = ncread(image_in_file, 'mask_manual');
mask_manual = logical(mask_manual);

ini = ncread(image_in_file, 'img', [1, 1, ini_index],   [inf, inf, 1]); 
ini = double(ini);

fin = ncread(image_in_file, 'img', [1, 1, ini_index+gap], [inf, inf, 1]); 
fin = double(fin);

ini_mask = ncread(image_in_file, 'mask_auto', [1, 1, ini_index],   [inf, inf, 1]); 
ini_mask = logical(ini_mask) & mask_manual;

fin_mask = ncread(image_in_file, 'mask_auto', [1, 1, ini_index+gap],   [inf, inf, 1]); 
fin_mask = logical(fin_mask) & mask_manual;

xw = ncread(image_in_file, 'x');
xw = double(xw);

yw = ncread(image_in_file, 'y');
yw = double(yw);

%% Display image and mask data

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

%% Set parameters, run PIV, and compute strain

% set parameters
samp_len = 60;
samp_spc = 30;
intr_len = 90;
num_pass = 1;
valid_radius = 120;
valid_max = 2;
valid_eps = 0.01;
spline_tension = 0.95;
min_frac_data = 0.50; 
min_frac_overlap = 0.20;

piv_result = piv(ini, fin, ini_mask, fin_mask, xw, yw, samp_len, samp_spc, ...
        intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
        spline_tension, min_frac_data, min_frac_overlap, true);         

% extract key variables from the results struct
xx = piv_result.x_grd(1,:); 
yy = piv_result.y_grd(:,1);
uu = piv_result.u_grd;
vv = piv_result.v_grd;
roi = piv_result.roi_grd;    
    
% compute strain
strain_result = post_strain(xx, yy, uu, vv, roi, 'nearest');

fprintf('DONE\n');

%% Display results

close all

% get data limits
mask_uv = ~isnan(uu) & ~isnan(vv);
[xxg, yyg] = meshgrid(xx, yy);
uv_xlim = [min(xxg(mask_uv)), max(xxg(mask_uv))];
uv_ylim = [min(yyg(mask_uv)), max(yyg(mask_uv))];

% plot displacements
%... displacement magnitude and direction, [mm/step]
figure
subplot(3,1,1);
mm = sqrt(uu.^2 + vv.^2);
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

% plot strain
%... Dd
figure
subplot(3,1,1);
imagesc(xx, yy, strain_result.Dd, 'AlphaData', ~isnan(strain_result.Dd)); 
colorbar;
set(gca, 'YDir', 'normal', 'XLim', uv_xlim, 'YLim', uv_ylim);
title('Dd');
%... Dv
subplot(3,1,2);
imagesc(xx, yy, strain_result.Dv, 'AlphaData', ~isnan(strain_result.Dv)); 
colormap(gca, flipud(colormap)); % flow is in negative x direction
colorbar;
set(gca, 'YDir', 'normal', 'XLim', uv_xlim, 'YLim', uv_ylim);
title('Dv');
%... spin
subplot(3,1,3);
imagesc(xx, yy, strain_result.spin, 'AlphaData', ~isnan(strain_result.spin)); 
colorbar;
set(gca, 'YDir', 'normal', 'XLim', uv_xlim, 'YLim', uv_ylim);
title('spin');

%% Analyze short series using piv_series()

junk_output_file = 'delete_me.nc';
step_range = [300, 303];

piv_series(junk_output_file, image_in_file, step_range, gap, samp_len, ...
    samp_spc, intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
    spline_tension, min_frac_data, min_frac_overlap, true)

%% Save parameters to mat file

save(param_out_file, 'samp_len',  'samp_spc', 'intr_len', 'num_pass', ...
    'valid_radius', 'valid_max', 'valid_eps', 'spline_tension', ...
    'min_frac_data', 'min_frac_overlap', 'gap');