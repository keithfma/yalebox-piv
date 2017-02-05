% Script. Scratch space for debugging the post-processing code suite

% parameters
exp_name = 'fault_ss_01_sidef';
data_dir = 'C:\Users\keithfma\home\dissertation\yalebox-exp-fault\data';
piv_netcdf = fullfile(data_dir, sprintf('%s_piv.nc', exp_name));
post_netcdf = fullfile(data_dir, sprintf('%s_post.nc', exp_name));
pro_bbox = [0.56, 0.0, 0.05, 0.1]; % limited in x, whole section in y
retro_bbox = [-0.5, 0.0, 0.1, 0.1];

% cleanup previous run
delete(post_netcdf);

% try series post-processing 
post_series(piv_netcdf, post_netcdf, pro_bbox, retro_bbox);
