% Script. Scratch space for debugging the post-processing code suite

% parameters
data_dir = 'C:\Users\keithfma\home\dissertation\yalebox-exp-fault\data';
piv_netcdf = fullfile(data_dir, 'fault_ss_01_siden_piv.nc');
post_netcdf = 'junk.nc';

% cleanup previous run
delete(post_netcdf);

% try series post-processing 
post_series(piv_netcdf, post_netcdf);
