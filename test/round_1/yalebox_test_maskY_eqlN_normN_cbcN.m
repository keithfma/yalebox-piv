% Script. Test case, round 1.
%
% Masked = Y
% Equalized histogram = N
% Use normalized xcorr = N
% Use CBC = N
         
% define parameters
ini_file = {'fault_ss_01_sidef_030.png', 'fault_ss_01_sidef_250.png'}; 
ini_mask_file = {'030_mask.mat', '250_mask.mat'};
fin_file = {'fault_ss_01_sidef_031.png', 'fault_ss_01_sidef_251.png'};
fin_mask_file = {'031_mask.mat', '251_mask.mat'};
out_file = {'030_031_maskY_eqlN_normN_cbcN.mat', '250_251_maskY_eqlN_normN_cbcN.mat'};
coords_file = 'coords.mat';
npass = 1;
samplen = 30;
sampspc = 30;
umax = 0.01; 
umin = -0.02;
vmax = 0.01;
vmin = -0.01;
ncbc = 1;
verbose = 1;
use_normxcorr2 = 0;

% load coordinates
load(coords_file, 'x', 'y');

for i = 1:2
    
    % read images and apply masks
    im = rgb2hsv(imread(ini_file{i}));
    ini = im(:,:,3);
    load(ini_mask_file{i}, 'mask_manual', 'mask_auto');
    ini = ini.*mask_manual.*mask_auto;
    
    im = rgb2hsv(imread(fin_file{i}));
    fin = im(:,:,3);
    load(fin_mask_file{i}, 'mask_manual', 'mask_auto');
    fin = fin.*mask_manual.*mask_auto;
    
    % % add a tiny bit of noise, to avoid a "template cannot all be the same" error
    % ini = ini-1e-6*rand(size(ini));
    % fin = fin-1e-6*rand(size(fin));    
    
    % run piv
    [xx, yy, uu, vv] = yalebox_piv_step(...
        ini, fin, x, y, npass, samplen, sampspc, umax, umin, vmax, vmin, ...
        ncbc, verbose, use_normxcorr2);
    
    [xxgrid, yygrid] = meshgrid(xx, yy);
    
    [displacement, spin, Dv, Dd, D2x, D2y, WkStar, AkStar] = ...
        yalebox_decompose_step(xxgrid, yygrid, uu, vv, ~isnan(uu));
    
    % save results to file
    save(out_file{i}, 'xx','yy', 'uu', 'vv', 'npass', 'samplen', 'sampspc', ...
        'umax', 'umin', 'vmax', 'vmin', 'ncbc', 'use_normxcorr2', ...
        'displacement', 'spin', 'Dv', 'Dd', 'D2x', 'D2y', 'WkStar', 'AkStar');
end