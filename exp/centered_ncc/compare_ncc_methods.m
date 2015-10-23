% Script to explore Mark Brandon's idea of setting masked regions equal to the
% window mean prior to computing the normalized cross correlation. As a
% shorthand, I call this "centering". The goal is to reduce the influence of
% masked regions on the fit. I test two verions: a crude version in which the
% interogation window is centered only once, and (2) a complete version in which
% the subimage for each shift is centered. Results are compared with the
% standard normalized cross correlation, and the masked normalized cross
% correlation.
%
% Cleaned-up for sharing.

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

% masked ncc, weighted
nxcr_masked_weighted = nxcr_masked.*(noverlap/max(noverlap(:)));
warning on images:removing:function

% centered ncc, simple 
nxcr_centered_simple = center_simple_normxcorr2(samp, intr);

% centered ncc, complete
nxcr_centered = center_normxcorr2(samp, intr);


%% Compare results

% define parameters
clim = [-1, 1];
pos = [0.05, 0.05, 0.9, 0.9];

% plot
figure('units', 'normalized', 'position', pos);

subplot(3,2,1)
imagesc(xcr);
title('Mathworks CC');
%caxis(clim);
colorbar
axis equal
axis tight

subplot(3,2,2)
imagesc(nxcr);
title('Mathworks NCC');
colorbar
caxis(clim);
axis equal
axis tight

subplot(3,2,3)
imagesc(nxcr_masked);
title('Padfield 2010 Masked NCC')
colorbar
caxis(clim);
axis equal
axis tight

subplot(3,2,4)
imagesc(nxcr_masked_weighted);
title('Padfield 2010 Masked NCC, Weighted')
colorbar
caxis(clim);
axis equal
axis tight

subplot(3,2,5)
imagesc(nxcr_centered_simple);
title('Centered NCC, Simple')
colorbar
caxis(clim);
axis equal
axis tight

subplot(3,2,6)
imagesc(nxcr_centered);
title('Centered NCC, Complete')
colorbar
caxis(clim);
axis equal
axis tight
