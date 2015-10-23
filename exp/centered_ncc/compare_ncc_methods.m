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
tic;
xcr = xcorr2(samp, intr);
t_xcr = toc;

% standard ncc
tic;
nxcr = normxcorr2(samp, intr);
t_nxcr = toc;

% masked ncc
warning off images:removing:function
tic;
[nxcr_masked, noverlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);
t_nxcr_masked = toc;
warning on images:removing:function

% masked ncc, weighted
tic;
nxcr_masked_weighted = nxcr_masked.*(noverlap/max(noverlap(:)));
t_nxcr_masked_weighted = toc+t_nxcr_masked;

% centered ncc, simple 
tic;
nxcr_centered_simple = center_simple_normxcorr2(samp, intr);
t_nxcr_centered_simple= toc;

% centered ncc, complete
tic;
nxcr_centered = center_normxcorr2(samp, intr);
t_nxcr_centered = toc;


%% Compare results

% define parameters
clim = [-1, 1];
pos = [0.05, 0.05, 0.9, 0.9];

% plot
figure('units', 'normalized', 'position', pos);

subplot(3,2,1)
imagesc(xcr);
title(sprintf('Mathworks CC, %.0f ms', 1000*t_xcr));
%caxis(clim);
colorbar
axis equal
axis tight

subplot(3,2,2)
imagesc(nxcr);
title(sprintf('Mathworks NCC, %.0f ms', 1000*t_nxcr));
colorbar
caxis(clim);
axis equal
axis tight

subplot(3,2,3)
imagesc(nxcr_masked);
title(sprintf('Padfield 2010 Masked NCC, %.0f ms', 1000*t_nxcr_masked));
colorbar
caxis(clim);
axis equal
axis tight

subplot(3,2,4)
imagesc(nxcr_masked_weighted);
title(sprintf('Padfield 2010 Masked NCC, Weighted, %.0f ms', 1000*t_nxcr_masked_weighted));
colorbar
caxis(clim);
axis equal
axis tight

subplot(3,2,5)
imagesc(nxcr_centered_simple);
title(sprintf('Centered NCC, Simple, %.0f ms', 1000*t_nxcr_centered_simple));
colorbar
caxis(clim);
axis equal
axis tight

subplot(3,2,6)
imagesc(nxcr_centered);
title(sprintf('Centered NCC, Complete, %.0f ms', 1000*t_nxcr_centered));
colorbar
caxis(clim);
axis equal
axis tight
