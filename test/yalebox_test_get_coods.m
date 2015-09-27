% Script. Create coordinate vectors for the two test datasets

% define parameters
wocofile = 'fault_ss_01_woco_sidef.png';
npts = 14;
outfile = 'coords.mat';

% setup environment
addpath('../');

% get coordinates and save results to file
[x, y, scale, offset] = yalebox_prep_world_coord(wocofile, npts, true);
save(outfile, 'x', 'y', 'scale', 'offset');
