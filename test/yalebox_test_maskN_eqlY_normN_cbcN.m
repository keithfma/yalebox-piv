% Script. Test case, round 1.
%
% Masked = N
% Equalized histogram = Y
% Use normalized xcorr = N
% Use CBC = N
         
% define parameters
ini_file = {'030_eql.mat', '250_eql.mat'}; 
fin_file = {'031_eql.mat', '251_eql.mat'};
out_file = {'030_031_maskN_eqlY_normN_cbcN.mat', '250_251_maskN_eqlY_normN_cbcN.mat'};
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


% read images, run PIV, and save results, x2
for i = 1:2
    
    % read images
    load(ini_file{i}, 'eql');
    ini = eql;
    load(fin_file{i}, 'eql');
    fin = eql;
    
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