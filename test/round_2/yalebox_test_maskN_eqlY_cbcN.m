% Script. Test case, round 2
%
% Masked = N
% Equalized histogram = Y
% Use CBC = N
         
% define parameters
ini_file = {'../030_eql.mat', '../250_eql.mat'}; 
fin_file = {'../031_eql.mat', '../251_eql.mat'};
ini_mask_file = {'../030_mask.mat', '../250_mask.mat'};
fin_mask_file = {'../031_mask.mat', '../251_mask.mat'};
out_file = {'030_031_maskN_eqlY_cbcN.mat', '250_251_maskN_eqlY_cbcN.mat'};
out_fig = {'uv_maskN_eqlY_cbcN.fig', 'dd_maskN_eqlY_cbcN.fig'};
fig_name = {'displacement: maskN_eqlY_cbcN', 'strain: maskN_eqlY_cbcN'};
coords_file = '../coords.mat';
npass = 1;
samplen = 30;
sampspc = 15;
umax = 0.015; 
umin = -0.025;
vmax = 0.015;
vmin = -0.015;
ncbc = 1;
verbose = 1;
use_normxcorr2 = 1;

% load coordinates
load(coords_file, 'x', 'y');

% read images, run PIV, and save results, x2
for i = 1:2
    
    % read images
    load(ini_file{i}, 'eql');
    ini = eql;
    load(fin_file{i}, 'eql');
    fin = eql;
    
    % add a tiny bit of noise, to avoid a "template cannot all be the same" error
    ini = min(1, max(0, ini-1e-6*rand(size(ini))));
    fin = min(1, max(0, fin-1e-6*rand(size(fin))));
    
    % read masks
    load(ini_mask_file{i}, 'mask_manual', 'mask_auto');
    ini_mask = mask_manual & mask_auto;
    load(fin_mask_file{i}, 'mask_manual', 'mask_auto');
    fin_mask = mask_manual & mask_auto;       
    
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
    
    % interpolate mask to output resolution
    [xgrid, ygrid] = meshgrid(x, y);
    mm = interp2(xgrid, ygrid, double(ini_mask & fin_mask), xxgrid, yygrid) > 0.5;

    % generate standard results image
    h1 = figure(1);    
    set(gcf, 'Name', fig_name{1});
    subplot(2, 1, i)
    imagesc(xx, yy, displacement.*mm);
    caxis([0, 0.01])
    colorbar
    axis equal
    axis tight
    hold on
    quiver(xx, yy, uu, vv, 4, '-k');
    hold off
    set(gca, 'YDir', 'normal')
    
    h2 = figure(2);
    set(gcf, 'Name', fig_name{2});
    subplot(2, 1, i)
    title('strain (Dd)')
    imagesc(xx, yy, Dd.*mm)
    caxis([0, 0.08])
    colorbar
    axis equal
    axis tight
    set(gca, 'YDir', 'normal')
    
end

% save figures
savefig(h1, out_fig{1});
savefig(h2, out_fig{2});
