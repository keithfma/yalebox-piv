% Scratch-paper script

%% Explore padding PIV velocity fields prior to deformation analysis

% parameters
data_file = '../yalebox-exp-fault/data/fault_ss_01_siden_piv_noloess.nc';
step = 360.5;

% load data
s = double(ncread(data_file, 'step'));
step_idx = find(s==step);
uu = double(ncread(data_file, 'u', [1, 1, step_idx], [inf, inf, 1]));
vv = double(ncread(data_file, 'v', [1, 1, step_idx], [inf, inf, 1]));
roi = logical(ncread(data_file, 'roi', [1, 1, step_idx], [inf, inf, 1]));
x = double(ncread(data_file, 'x'));
y = double(ncread(data_file, 'y'));

% compute instantaneous deformation (includes padding)
[L, F, S] = post_strain(x, y, uu, vv, roi, 'nearest');