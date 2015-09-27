% Script. Test case, round 3
%
% Equalized histogram = N
% Use CBC = N
         
% define parameters
ini_file = {'../030_eql.mat', '../250_eql.mat'}; 
fin_file = {'../031_eql.mat', '../251_eql.mat'};
ini_mask_file = {'../030_mask.mat', '../250_mask.mat'};
fin_mask_file = {'../031_mask.mat', '../251_mask.mat'};
out_file = {'030_031_eqlY_cbcN.mat', '250_251_eqlY_cbcN.mat'};
out_fig = {'uv_eqlY_cbcN.fig', 'dd_eqlY_cbcN.fig'};
fig_name = {'displacement: eqlY_cbcN', 'strain: eqlY_cbcN'};
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
    
    % read masks
    load(ini_mask_file{i}, 'mask_manual', 'mask_auto');
    ini_mask = mask_manual & mask_auto;
    load(fin_mask_file{i}, 'mask_manual', 'mask_auto');
    fin_mask = mask_manual & mask_auto;    
    
    % apply masks
    ini = ini.*ini_mask;
    fin = fin.*fin_mask;
    
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

    % generate standard results image
    h1 = figure(1);    
    set(gcf, 'Name', fig_name{1});
    subplot(2, 1, i)
    imagesc(xx, yy, displacement);
    caxis([0, 0.01])
    colorbar
    axis equal
    axis tight
    hold on
    quiver(xx(1:2:end), yy(1:2:end), ...
        uu(1:2:end, 1:2:end), vv(1:2:end, 1:2:end),...
        4, '-k');
    hold off
    set(gca, 'YDir', 'normal')
    
    h2 = figure(2);
    set(gcf, 'Name', fig_name{2});
    subplot(2, 1, i)
    title('strain (Dd)')
    imagesc(xx, yy, Dd)
    caxis([0, 0.06])
    colorbar
    axis equal
    axis tight
    set(gca, 'YDir', 'normal')
    
end

% save figures
savefig(h1, out_fig{1});
savefig(h2, out_fig{2});