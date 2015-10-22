% Script to explore Mark Brandon's idea of setting masked regions equal to the
% window mean prior to computing the normalized cross correlation. As a
% shorthand, I call this "centering". The goal is to reduce the influence of
% masked regions on the fit. I test two verions: a crude version in which the
% interogation window is centered only once, and (2) a complete version in which
% the subimage for each shift is centered. Results are compared with the
% standard normalized cross correlation, and the masked normalized cross
% correlation.

%% Initialize

% define parameters

% data_file = 'simple_shear_interior.mat';
data_file = 'simple_shear_lower_bnd.mat';

% load data
load(data_file, 'samp', 'intr');

%% Run cases

% standard cc
xcr = xcorr2(samp, intr);

% standard ncc
nxcr = normxcorr2(samp, intr);

% masked ncc
warning off images:removing:function
[nxcr_masked, noverlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);

% masked ncc, trimmed
nxcr_masked_trimmed = nxcr_masked;
nxcr_masked_trimmed(noverlap<0.75*max(noverlap(:))) = 0;
warning on images:removing:function

% centered ncc
nxcr_centered = center_simple_normxcorr2(samp, intr);

%% Compare results

% define parameters
clim = [-1, 1];
pos = [0.05, 0.2, 0.9, 0.5];

% plot
figure('units', 'normalized', 'position', pos);

subplot(1,5,1)
imagesc(xcr);
title('Stardard CC');
%caxis(clim);
colorbar
axis equal
axis tight

subplot(1,5,2)
imagesc(nxcr);
title('Stardard NCC');
colorbar
caxis(clim);
axis equal
axis tight

subplot(1,5,3)
imagesc(nxcr_masked);
title('Masked NCC')
colorbar
caxis(clim);
axis equal
axis tight

subplot(1,5,4)
imagesc(nxcr_masked_trimmed);
title('Masked NCC, Trimmed')
colorbar
caxis(clim);
axis equal
axis tight

subplot(1,5,5)
imagesc(nxcr_centered);
title('Centered NCC')
colorbar
caxis(clim);
axis equal
axis tight