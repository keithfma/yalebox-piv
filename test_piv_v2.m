function [xx, yy, uu_piv, vv_piv, uu_exact, vv_exact] = test_piv_v2(varargin)
% function [] = test_piv(varargin)
%
% Test PIV performance using case using synthetic image pair
%
% NOTE: This is a "high-level" test that covers all of the PIV subsystems.
%
% Optional Arguments ('Name', Value):
%   'image_file': String, name of pre-processed image file, as produced by
%       prep_series(), to read for raw image, default = ./test/default_image.nc
%   'image_index': Integer, (1-based) index of image in image file to use for
%       raw image, default = 1
%   'image_pos': 4-element position vector indicating the limits of the image to
%       extract and deform. Must contain only sand (all within the ROI),
%       in meters, default = [-0.12, 0.005, 0.092, 0.07]
%   'u1': 2-element vector specifying first end-member translation in pixels,
%       default = 20*[cosd(45), -sind(45)]
%   'u2': 2-element vector specifying second end-member translation in pixels,
%       default = -20*[cosd(45), -sind(45)]
%   'theta': Scalar, orientation of shear band specified as
%       counter-clockwise angle to the positive x-axis, in degrees, limited to
%       range 0 - 90, default = 45
%   'shear_width': Scalar, width of shear band, as fraction of image width,
%       default = 0.25   
%   'pad_width': Scalar, width of edge padding to add to image (to accomodate
%       edge displacements) as a fraction of image size, default = 0.1
%   'samp_len': piv() parameter, default [60, 30]
%   'samp_spc': piv() parameter, default 15
%   'intr_len': piv() parameter, default [90, 40]
%   'num_pass': piv() parameter, default [1, 2]
%   'valid_radius: piv() parameter, default 45
%   'valid_max': piv() parameter, default 2
%   'valid_eps': piv() parameter, default 0.1
%   'min_frac_data': piv() parameter, default 0.5
%   'min_frac_overlap': piv() parameter, default 0.25
%   'verbose': Scalar logical, set true to enable verbose reporting for all
%       components of the analysis
% %

%% parse arguments

% constants
src_dir = fileparts(mfilename('fullpath'));

ip = inputParser();

ip.addParameter('image_file', fullfile(src_dir, 'test', 'default_image.nc'), ...
    @(x) exist(x, 'file') == 2);
ip.addParameter('image_index', 1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
ip.addParameter('image_pos', [-0.12, 0.005, 0.092, 0.07], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 4}));
ip.addParameter('velocity_a',  20*[cosd(45), -sind(45)], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('velocity_b', -20*[cosd(45), -sind(45)], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('theta', 45, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=' 90}));
ip.addParameter('defm_width', 0.25, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
ip.addParameter('pad_width', 0.1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0}));
ip.addParameter('samp_len', [60, 30]); % validation handled by PIV routines
ip.addParameter('samp_spc', 15);
ip.addParameter('intr_len', [90, 40]);
ip.addParameter('num_pass', [1, 2]);
ip.addParameter('valid_radius', 45);
ip.addParameter('valid_max', 2);
ip.addParameter('valid_eps', 0.1);
ip.addParameter('min_frac_data', 0.5);
ip.addParameter('min_frac_overlap', 0.25);
ip.addParameter('verbose', true, ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {'scalar', 'binary'}));

ip.parse(varargin{:});
args = ip.Results;

if args.verbose
    fprintf('%s: arguments:\n', mfilename);
    disp(args)
end


%% read and crop raw image

if args.verbose
    fprintf('%s: read and crop raw image\n', mfilename);
end

xw = double(ncread(args.image_file, 'x'));
yw = double(ncread(args.image_file, 'y'));
img = double(ncread(args.image_file, 'img', [1, 1, args.image_index], [inf, inf, 1]));
mask_auto = ncread(args.image_file, 'mask_auto', [1, 1, args.image_index], [inf, inf, 1]);
mask_manu = ncread(args.image_file, 'mask_manual');
roi = mask_auto & mask_manu;

min_col = find(xw >= args.image_pos(1), 1, 'first');
max_col = find(xw <= args.image_pos(1) + args.image_pos(3), 1, 'last');
min_row = find(yw >= args.image_pos(2), 1, 'first');
max_row = find(yw <= args.image_pos(2) + args.image_pos(4), 1, 'last');

img = img(min_row:max_row, min_col:max_col);
roi = roi(min_row:max_row, min_col:max_col);


if any(~roi(:))
    error('%s: image_pos limits must include only sand (ROI)', mfilename);
end

% pad image boundaries (and coordinates) to accomodate edge displacements
if args.verbose
    fprintf('%s: pad image to accomodate edge displacements\n', mfilename);
end

pad_dim = ceil(args.pad_width*size(img));
img = padarray(img, pad_dim, 0, 'both');
roi = padarray(roi, pad_dim, 0, 'both');

% create new coordinate vectors, pixel units with origin at padded image center
[x_img, y_img] = meshgrid(1:size(img, 2), 1:size(img, 1));
x_img = x_img - mean(x_img(:));
y_img = y_img - mean(y_img(:));

%% generate synthetic images

% compute exact displacement field for specified displacements and boundary
if args.verbose
    fprintf('%s: compute exact displacement field\n', mfilename);
end

[u_img, v_img] = exact_velocity(...
    x_img, y_img, args.defm_width/2*size(img, 2), args.theta, ...
    args.velocity_a, args.velocity_b);

if args.verbose
    fprintf('%s: generate synthetic images\n', mfilename);
end

% deform images
ini = imwarp(img, 0.5*cat(3, u_img, v_img), 'cubic'); % dir is ok
fin = imwarp(img, -0.5*cat(3, u_img, v_img), 'cubic');

% deform masks
ini_roi = imwarp(double(roi), 0.5*cat(3, u_img, v_img), 'cubic');
ini_roi = logical(round(ini_roi));
fin_roi = imwarp(double(roi), -0.5*cat(3, u_img, v_img), 'cubic');
fin_roi = logical(round(fin_roi));

% enforce limits by threshold
% NOTE: could stretch limits instead, unclear if this matters
% NOTE: add a tiny offset so that no sand pixels are exactly zero
tiny = 1e-5;
ini(ini < 0) = tiny;
ini(ini > 1) = 1;
fin(fin < 0) = tiny;
fin(fin > 1) = 1;

% reapply mask
ini(~ini_roi) = 0;
fin(~fin_roi) = 0;

% display initial and final synthetic images
figure('Position', get(0, 'ScreenSize'))

subplot(1,2,1)
imagesc([x_img(1), x_img(end)], [y_img(1), y_img(end)], ini);
set(gca, 'YDir', 'normal', 'XGrid', 'on', 'YGrid', 'on', 'GridColor', 'w');
axis equal tight
title('Initial Synthetic Image')

subplot(1,2,2)
imagesc([x_img(1), x_img(end)], [y_img(1), y_img(end)], fin);
set(gca, 'YDir', 'normal', 'XGrid', 'on', 'YGrid', 'on', 'GridColor', 'w');
axis equal tight
title('Final Synthetic Image')

%% run PIV analysis on synthetic images, get exact solution on same grid

if args.verbose
    fprintf('%s: run PIV analysis\n', mfilename);
end

result = piv(...     
    ini, fin, ini_roi, fin_roi, x_img(1,:), y_img(:,1), args.samp_len, args.samp_spc, ...
    args.intr_len, args.num_pass, args.valid_radius, args.valid_max, ...
    args.valid_eps, args.min_frac_data, args.min_frac_overlap, ...
    args.verbose);

% repackage key results
coord_obs = struct('x', result.x_grd, 'y', result.y_grd, 'roi', result.roi_grd);
velocity_obs = struct('u', result.u_grd, 'v', result.v_grd, ...
    'm', sqrt(result.u_grd.^2 + result.v_grd.^2), ...
    'theta', atan2d(result.v_grd, result.u_grd));

% apply roi mask to velocity fields
velocity_obs.u(~coord_obs.roi) = NaN;
velocity_obs.v(~coord_obs.roi) = NaN;

% compute exact velocities and package results
[u_ext, v_ext] = exact_velocity(...
    coord_obs.x, coord_obs.y, args.defm_width/2*size(img, 2), args.theta, ...
    args.velocity_a, args.velocity_b);
velocity_ext = struct('u', u_ext, 'v', v_ext, ...
    'm', sqrt(u_ext.^2 + v_ext.^2), ...
    'theta', atan2d(v_ext, u_ext));

%% compute velocity errors and store in struct

velocity_error = struct();
velocity_fields = fieldnames(velocity_ext);
for ii = 1:numel(velocity_fields)
    fn = velocity_fields{ii};
    velocity_error.(fn) = velocity_ext.(fn) - velocity_obs.(fn);
end
velocity_error_abs = structfun(@abs, velocity_error, 'UniformOutput', false);

%% run deformation analysis, get exact deformation parameters on same grid

strain_obs = post_strain(coord_obs.x(1,:), coord_obs.y(:,1), velocity_obs.u, ...
                         velocity_obs.v, coord_obs.roi, 'nearest');

strain_ext = exact_strain(...
    coord_obs.x, coord_obs.y, args.defm_width/2*size(img, 2), args.theta, ...
    args.velocity_a, args.velocity_b);

%% compute strain errors and store in struct

strain_error = struct();
strain_fields = fieldnames(strain_ext);
for ii = 1:numel(strain_fields)
    fn = strain_fields{ii};
    strain_error.(fn) = strain_ext.(fn) - strain_obs.(fn);
end
strain_error_abs = structfun(@abs, strain_error, 'UniformOutput', false);

%% IMPORTANT

% It looks like there is a sign error somewhere spin observed and exact are off
% by factor of -1

%% START DEBUG

keyboard

%% try clustering to split error variables

% TODO: explore k-medoid clustering with Mahalnobis distance, if needed

num_clusters = 3;
num_replicates = 25;

% DEBUG: add coordinates to clustering
velocity_error.x = coord_obs.x;
velocity_error.y = coord_obs.y;
velocity_error.d = coord_obs.x*sind(args.theta) + coord_obs.y*cosd(args.theta);

% DEBUG: add distance to shear zone center to clustering


velocity_cluster = cluster_errors(velocity_error, {'u', 'v', 'd'}, ...
    coord_obs.roi, num_clusters, num_replicates, true);
 
% % cluster select velocity fields
% velocity_cluster_fields = {'u', 'v', 'm'};
% velocity_features = nan(sum(coord_obs.roi(:)), numel(velocity_cluster_fields));
% for ii = 1:numel(velocity_cluster_fields)
%     fn = velocity_cluster_fields{ii};
%     velocity_features(:,ii) = velocity_error_abs.(fn)(coord_obs.roi);
% end
% velocity_cluster = nan(size(coord_obs.roi));
% velocity_cluster(coord_obs.roi) = kmeans(velocity_features, num_clusters, ...
%     'Replicates', num_replicates);

subplot(1,2,1)
imagesc(velocity_error.m);
subplot(1,2,2)
imagesc(velocity_cluster);


function clusters = cluster_errors(err, fields, roi, num_clusters, num_replicates, whiten)
%
% Cluster error data using all specified fields as features
%
% Arguments:
%   err: Struct, error data to be clustered, only specified fields are used
%   fields: Cell vector, fieldnames to include as features, set to
%       fieldnames(err) to use all available fields
%   roi: 2D matrix, logical array indicating data footprint
%   num_clusters: Scalar integer, number of clusters to fit
%   num_replicates: Scalar integer, number of replicate clustering attempts,
%       helps avoid local minima
%   whiten: Scalar logical, set true to apply PCA whitening transform
%
% Returns:
%   2D matrix of cluster ids at each point
%
% NOTE: PCA whitening transforms data such that they have covariance matrix ==
%   identity matrix, meaning variable variances are 1 and variables do not
%   covary. This is analagous to the Mahalnobis distance metric. Sounds good,
%   but does not always improve performance.
% %

% build feature matrix
features = nan(sum(roi(:)), numel(fields));
for ii = 1:size(features, 2)
    features(:,ii) = err.(fields{ii})(roi);
end

% apply PCA whitening transform
if whiten
    [~, score, latent] = pca(features);
    features = score./repmat(sqrt(latent)' + 0.1, size(score, 1), 1);
end

% find clusters and store in correct positions
clusters = nan(size(roi));
clusters(roi) = kmeans(features, num_clusters, 'Replicates', num_replicates);

return

% % cluster all error fields, add results to the struct
% fields = fieldnames(err); % err_abs has same fields
% for ii = 1:numel(fields)
%     fn = fields{ii};
    

% % get cluster IDs
% num_clusters = 2;
% roi = ~isnan(abs_err.Dd); % is this the same for all variables?
% var = abs_err.Dd(roi);
% cluster_id = nan(size(roi));
% cluster_id(roi) = kmeans(var, num_clusters);

% % plot each cluster and its histogram
% for cc = 1:num_clusters
%     figure;
%     
%     rng = [min(
%     subplot(2, 1, 1);
%     imagesc(


% %% display exact and measured displacement fields
% 
% figure('Position', get(0, 'ScreenSize'))
% 
% num_vec = 25; % desired num quiver vectors along largest dim, for downsampling
% dfact = floor(max(size(x_piv))/num_vec);
% 
% subplot(1, 2, 1)
% m_piv_ext = sqrt(u_piv_ext.^2 + v_piv_ext.^2);
% imagesc(x_piv(1,:), y_piv(:,1), m_piv_ext, 'AlphaData', roi_piv_obs);
% set(gca, 'YDir', 'Normal', 'XGrid', 'on', 'YGrid', 'on', 'GridColor', 'w');
% hold on;
% quiver(x_piv(1:dfact:end, 1:dfact:end), ...
%     y_piv(1:dfact:end, 1:dfact:end), ...
%     u_piv_ext(1:dfact:end, 1:dfact:end), ...
%     v_piv_ext(1:dfact:end, 1:dfact:end), '-k');
% cb = colorbar;
% cb.Label.String = 'Displacement Magnitude [pixels]';
% axis equal tight
% title('Exact Displacement Magnitude and Direction')
% 
% subplot(1, 2, 2)
% m_piv_obs = sqrt(u_piv_obs.^2 + v_piv_obs.^2);
% imagesc(x_piv(1,:), y_piv(:,1), m_piv_obs, 'AlphaData', roi_piv_obs);
% set(gca, 'YDir', 'Normal', 'XGrid', 'on', 'YGrid', 'on', 'GridColor', 'w');
% hold on;
% quiver(x_piv(1:dfact:end, 1:dfact:end), ...
%     y_piv(1:dfact:end, 1:dfact:end), ...
%     u_piv_obs(1:dfact:end, 1:dfact:end), ...
%     v_piv_obs(1:dfact:end, 1:dfact:end), '-k');
% cb = colorbar;
% cb.Label.String = 'Displacement Magnitude [pixels]';
% axis equal tight
% title('PIV Displacement Magnitude and Direction')
% 
% %% Analyze veocity errors
% 
% if args.verbose
%     fprintf('%s: error analysis\n', mfilename);
% end
% 
% u_error = u_piv_ext - u_piv_obs;
% v_error = v_piv_ext - v_piv_obs;
% m_error = sqrt(u_error.^2 + v_error.^2);
% theta_error = atan2d(v_piv_ext, u_piv_ext) - atan2d(v_piv_obs, u_piv_obs);
% 
% % print error quantiles
% qnt = 0 : 0.10 : 1;
% u_abs_error_qnt = quantile(abs(u_error(:)), qnt);
% v_abs_error_qnt = quantile(abs(v_error(:)), qnt);
% m_error_qnt = quantile(sqrt(u_error(:).^2 + v_error(:).^2), qnt);
% 
% fprintf('\n---- \n');
% fprintf('Displacement Vector Error Quantiles\n\n');
% fprintf('Quantile\t| U Abs Error\t| V Abs Error\t| Mag Error\n'); 
% for i = 1:length(qnt)
%     fprintf('%.2f\t| %.2e\t| %.2e\t| %.2e\n', ...
%         qnt(i), u_abs_error_qnt(i), v_abs_error_qnt(i), m_error_qnt(i));
% end
% fprintf('---- \n\n');
% 
% % plot error maps
% figure('Position', get(0, 'ScreenSize'));
% 
% subplot(2, 2, 1);
% imagesc(x_piv(1,:), y_piv(:,1), u_error, 'AlphaData', roi_piv_obs);
% set(gca, 'YDir', 'normal');
% axis equal tight
% cb = colorbar;
% cb.Label.String = 'Error [pixels]';
% title('U Error')
% 
% subplot(2, 2, 2);
% imagesc(x_piv(1,:), y_piv(:,1), v_error, 'AlphaData', roi_piv_obs);
% set(gca, 'YDir', 'normal');
% axis equal tight
% cb = colorbar;
% cb.Label.String = 'Error [pixels]';
% title('V Error')
% 
% subplot(2, 2, 3);
% imagesc(x_piv(1,:), y_piv(:,1), m_error, 'AlphaData', roi_piv_obs);
% set(gca, 'YDir', 'normal');
% axis equal tight
% cb = colorbar;
% cb.Label.String = 'Error [pixels]';
% title('Magnitude Error')
% 
% subplot(2, 2, 4);
% imagesc(x_piv(1,:), y_piv(:,1), theta_error, 'AlphaData', roi_piv_obs);
% set(gca, 'YDir', 'normal');
% axis equal tight
% cb = colorbar;
% cb.Label.String = 'Error [degrees]';
% title('Theta Error')
% 
% % plot error histograms
% num_bins = 50;
% 
% figure('Position', get(0, 'ScreenSize'))
% 
% subplot(2, 2, 1)
% hist(u_error(~isnan(u_error)), num_bins);
% set(gca, 'XLim', [min(u_error(:)), max(u_error(:))]);
% hh = findobj(gca, 'Type', 'patch');
% hh.LineStyle = 'none';
% title('U Error')
% xlabel('Error [pixels]');
% ylabel('Count [pixels]');
% 
% subplot(2, 2, 2)
% hist(v_error(~isnan(v_error)), num_bins);
% set(gca, 'XLim', [min(v_error(:)), max(v_error(:))]);
% hh = findobj(gca, 'Type', 'patch');
% hh.LineStyle = 'none';
% title('V Error')
% xlabel('Error [pixels]');
% ylabel('Count [pixels]');
% 
% subplot(2, 2, 3)
% hist(m_error(~isnan(m_error)), num_bins);
% set(gca, 'XLim', [min(m_error(:)), max(m_error(:))]);
% hh = findobj(gca, 'Type', 'patch');
% hh.LineStyle = 'none';
% title('Magnitude Error')
% xlabel('Error [pixels]');
% ylabel('Count [pixels]');
% 
% subplot(2, 2, 4)
% hist(theta_error(~isnan(theta_error)), num_bins);
% set(gca, 'XLim', [min(theta_error(:)), max(theta_error(:))]);
% hh = findobj(gca, 'Type', 'patch');
% hh.LineStyle = 'none';
% title('Theta Error')
% xlabel('Error [degrees]');
% ylabel('Count [pixels]');

%% DEBUG 

function [uu, vv] = exact_velocity(xx, yy, hh, theta, uv_a, uv_b)
% Compute analytical velocity at points in xx, yy
%
% Arguments:
%   xx, yy: Matrices, coordinates at which to compute velocity
%   hh: Scalar, half-width of deformation zone, units match xx, yy
%   theta: Scalar, rotation angle for the deformation zone
%   uv_a, uv_b: 2-element vectors, end-member velocity vectors
%
% Returns:
%   uu, vv: Matrices, x- and y-components of velocity at points in xx, yy
% %

yr = xx*sind(theta) + yy*cosd(theta);
ww = nan(size(xx));
ww(yr > hh) = 1;
ww(yr < -hh) = 0;
defm = yr >= -hh & yr <= hh;
ww(defm) = 0.5 + yr(defm)/(2*hh);
uu = uv_a(1)*ww + uv_b(1)*(1 - ww);
vv = uv_a(2)*ww + uv_b(2)*(1 - ww);

return


function [dudx, dudy, dvdx, dvdy] = exact_gradient(xx, yy, hh, theta, uv_a, uv_b)
% Compute analytical deformation quantities at points in xx, yy
%
% Arguments:
%   xx, yy: Matrices, coordinates at which to compute velocity
%   hh: Scalar, half-width of deformation zone, units match xx, yy
%   theta: Scalar, rotation angle for the deformation zone
%   uv_a, uv_b: 2-element vectors, end-member velocity vectors
%
% Returns:
%   dudx, dudy, dvdx, dvdy: Matrices, spatial velocity gradients, a.k.a.
%       elements of the velocity gradient tensor
% %

yr = xx*sind(theta) + yy*cosd(theta);
defm = yr >= -hh & yr <= hh;

dudx = zeros(size(xx));
dvdx = zeros(size(xx));
dudy = zeros(size(xx));
dvdy = zeros(size(xx));

dudx(defm) = (uv_a(1) - uv_b(1))*sind(theta)/(2*hh);
dvdx(defm) = (uv_a(2) - uv_b(2))*sind(theta)/(2*hh);
dudy(defm) = (uv_a(1) - uv_b(1))*cosd(theta)/(2*hh);
dvdy(defm) = (uv_a(2) - uv_b(2))*cosd(theta)/(2*hh);

return


function strain = exact_strain(xx, yy, hh, theta, uv_a, uv_b)
% Compute analytical strain parameters at points in xx, yy
%

% 
% Arguments:
%   xx, yy: Matrices, coordinates at which to compute velocity
%   hh: Scalar, half-width of deformation zone, units match xx, yy
%   theta: Scalar, rotation angle for the deformation zone
%   uv_a, uv_b: 2-element vectors, end-member velocity vectors
% 
% NOTE: Solution for spatial velocity gradients is analytical, calculation of
%   the derived strain quantities is numerical. It may be possible to compute
%   these values directly with a bit of algebraic leg-work, but it is not
%   obviously worth it to do so (yet).
% %

% compute spatial velocity gradients (analytical)
[dudx, dudy, dvdx, dvdy] = exact_gradient(xx, yy, hh, theta, uv_a, uv_b);

% compute deformation gradient tensor: F = dx/dX = d/dX( X+u ) = I + du/dX
F11 = 1 + dudx;
F12 = dudy;
F21 = dvdx;
F22 = 1 + dvdy;

% preallocate derived strain parameters
[ny, nx] = size(dudx);
S1 = nan(ny, nx);
S2 = nan(ny, nx);
S1x = nan(ny, nx);
S1y = nan(ny, nx);
S2x = nan(ny, nx);
S2y = nan(ny, nx);
spin = nan(ny, nx);

% calculate derived strain parameters
for kk = 1:(nx*ny)
    F = [F11(kk), F12(kk); F21(kk), F22(kk)];
    % polar decomposition, F = VR, B = V^2 = FF', and R = (V^-1)F,
    % %     where B is the Left Cauchy-Green tensor.
    B = F*F';
    % % eigen solution sorted with eigenvalues in descending order
    [T, lambda] = eig(B);
    [lambda, order] = sort(diag(lambda),'descend');
    S1(kk) = sqrt(lambda(1));
    S2(kk) = sqrt(lambda(2));
    T = T(:, order);
    % store unit vectors for maximum and minimum extension rate direction
    S1x(kk) = T(1,1);
    S1y(kk) = T(2,1);
    S2x(kk) = T(1,2);
    S2y(kk) = T(2,2);
    % rotation matrix and spin
    V_inverse = T*diag(lambda.^(-1/2))*T';
    R = V_inverse*F;
    spin(kk) = atan2(R(2,1), R(1,1));
end
% equivalent steady rate-of-deformation / stretching tensor components
%   NOTE: D = logm(S), which is trivial for S in principle form
D1 = log(S1);
D2 = log(S2);

% calculate invariants, assuming plane strain (D3=0)
Dt = sqrt(D1.^2 + D2.^2);
Dd = sqrt(((D1-D2).^2 + D1.^2 + D2.^2)/3);
Dv = sqrt(1/3)*(D1 + D2); 
 
% calculate kinematic numbers
Wk      = 2*spin./(sqrt(2)*Dt);
Wk_star = 2*spin./(sqrt(2)*Dd);
Ak      = Dv./(sqrt(2)*Dt);
Ak_star = Dv./(sqrt(2)*Dd);

% copy output values to struct
strain = struct('F11', F11, 'F12', F12, 'F21', F21, 'F22', F22, 'S1', S1, ...
    'S2', S2, 'S1x', S1x, 'S1y', S1y, 'S2x', S2x, 'S2y', S2y, 'spin', spin, ...
    'D1', D1, 'D2', D2, 'Dt', Dt, 'Dd', Dd, 'Dv', Dv, 'Wk', Wk, 'Ak', Ak, ...
    'Wk_star', Wk_star, 'Ak_star', Ak_star);

return
