% Script. Create coordinate vectors for the two test datasets

% define parameters
wocofile_1 = '1/fault_ss_01_woco_sidef.png';
npts_1 = 14;
outfile_1 = '1/coords.mat';

wocofile_2 = '2/fault_ss_01_woco_siden.png';
npts_2 = 14;
outfile_2 = '2/coords.mat';

% setup environment
addpath('../');

% get coordinates and save results to file
% [x, y, scale, offset] = yalebox_prep_world_coord(wocofile_1, npts_1, true);
% save(outfile_1, 'x', 'y', 'scale', 'offset');

pause
close all

[x, y, scale, offset] = yalebox_prep_world_coord(wocofile_2, npts_2, true);
save(outfile_2, 'x', 'y', 'scale', 'offset');

