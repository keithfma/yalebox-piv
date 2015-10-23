% Script to run a test case and plot the results

%% parameters

data_file = 'simple_shear_interior.mat';


%% run

load(data_file, 'xcr');
[rpeak, cpeak, ok] = peak_optim_fourier(xcr);