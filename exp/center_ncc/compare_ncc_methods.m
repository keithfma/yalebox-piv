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

xcr = normxcorr2(samp, intr);

warning off images:removing:function
[xcr_masked, noverlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);
% xcr_masked(noverlap<0.5*max(noverlap(:))) = 0;
warning on images:removing:function

xcr_center_simple = center_simple_normxcorr2(samp, intr);

