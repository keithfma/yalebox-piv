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

% standard ncc
xcr = normxcorr2(samp, intr);

% masked ncc
warning off images:removing:function
[xcr_masked, noverlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);

% masked ncc, trimmed
xcr_masked_trimmed = xcr_masked;
xcr_masked_trimmed(noverlap<0.75*max(noverlap(:))) = 0;
warning on images:removing:function

% centered ncc
xcr_centered = center_simple_normxcorr2(samp, intr);

%% Compare results

% define parameters
clim = [-1, 1];
pos = [0.05, 0.2, 0.9, 0.5];

% plot
figure('units', 'normalized', 'position', pos);

subplot(1,4,1)
imagesc(xcr);
title('Stardard NCC');
caxis(clim);
axis equal
axis tight

subplot(1,4,2)
imagesc(xcr_masked);
title('Masked NCC')
caxis(clim);
axis equal
axis tight

subplot(1,4,3)
imagesc(xcr_masked_trimmed);
title('Masked NCC, Trimmed')
caxis(clim);
axis equal
axis tight

subplot(1,4,4)
imagesc(xcr_centered);
title('Centered NCC')
caxis(clim);
axis equal
axis tight




