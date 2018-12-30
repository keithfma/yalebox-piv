% Script. Explore and save piv parameters for experiment. Template to be copied
% and modified for specific experiments.

%% define and save parameters

image_file = 'XXX_image.nc';
piv_file = 'XXX_piv.nc';
param_file = 'XXX_piv_param.mat';
gap = 1;
step_range = [NaN, NaN]; % saved, but not used here
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

save(param_file, 'image_file', 'piv_file', 'gap', 'step_range', ...
    'samp_len',  'samp_spc', 'intr_len', 'num_pass', 'valid_radius', ...
    'valid_max', 'valid_eps', 'spline_tension', 'min_frac_data', ...
    'min_frac_overlap');

%% run analysis

ini_step = 350;

    [xw, yw, ini, ini_mask, fin, fin_mask] = ...
        util_read_input_pair(image_file, ini_step, gap);

piv_result = piv(ini, fin, ini_mask, fin_mask, xw, yw, samp_len, samp_spc, ...
        intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
        spline_tension, min_frac_data, min_frac_overlap, true);         

strain_result = post_strain(piv_result.x_grd(1,:), piv_result.y_grd(:,1), ...
    piv_result.u_grd, piv_result.v_grd, piv_result.roi_grd, 'nearest');

%% display results

quick_plot_image(xw, yw, ini, ini_mask, 'Initial: ');
quick_plot_image(xw, yw, fin, fin_mask, 'Final: ');
quick_plot_piv(piv_result);
quick_plot_strain(strain_result);