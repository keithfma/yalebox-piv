function [] = piv_test(varargin)
% function [] = piv_test(varargin)
%
% Test PIV performance using case using synthetic image pair. This is a
% high-level test that covers all of the PIV subsystems.
%
% Optional Arguments ('Name', Value):
% 
%   'image_file': String, name of pre-processed image file, as produced by
%       prep_series(), to read for raw image,
%       default = ./data/prep_testcase_0.mat
% 
%   'image_index': Integer, (1-based) index of image in image file to use for
%       raw image, default = 1
% 
%   'image_pos': 4-element position vector indicating the limits of the image to
%       extract and deform. Must contain only sand (all within the ROI),
%       in meters, default = [-0.12, 0.005, 0.092, 0.07]
% 
%   'u1': 2-element vector specifying first end-member translation in pixels,
%       default = 20*[cosd(45), -sind(45)]
% 
%   'u2': 2-element vector specifying second end-member translation in pixels,
%       default = -20*[cosd(45), -sind(45)]
% 
%   'theta': Scalar, orientation of shear band specified as
%       counter-clockwise angle to the positive x-axis, in degrees, limited to
%       range 0 - 90, default = 45
% 
%   'shear_width': Scalar, width of shear band, as fraction of image width,
%       default = 0.25   
% 
%   'pad_width': Scalar, width of edge padding to add to image (to accomodate
%       edge displacements) as a fraction of image size, default = 0.1
% 
%   'samp_len': piv() parameter, default [60, 30]
% 
%   'samp_spc': piv() parameter, default 15
% 
%   'intr_len': piv() parameter, default [90, 40]
% 
%   'num_pass': piv() parameter, default [1, 2]
% 
%   'valid_radius: piv() parameter, default 45
% 
%   'valid_max': piv() parameter, default 2
% 
%   'valid_eps': piv() parameter, default 0.1
% 
%   'min_frac_data': piv() parameter, default 0.5
% 
%   'min_frac_overlap': piv() parameter, default 0.25
% %

update_path('piv', 'post', 'akde');

% parse arguments  ----------

% constants
src_dir = fileparts(mfilename('fullpath'));

ip = inputParser();

ip.addParameter('image_file', fullfile(src_dir, 'data', 'prep_testcase_0.mat'), ...
    @(x) exist(x, 'file') == 2);
ip.addParameter('image_index', 1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
ip.addParameter('image_pos', [-0.1, 0.01, 0.1, 0.07], ...
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

ip.parse(varargin{:});
args = ip.Results;

fprintf('%s: arguments:\n', mfilename);
disp(args)


% read and crop raw image  ----------

fprintf('%s: read and crop raw image\n', mfilename);

image_data = matfile(args.image_file, 'Writable', false);
xw = image_data.x;
yw = image_data.y;
img = double(image_data.img(:, :, args.image_index));
mask_auto = image_data.mask_auto(:, :, args.image_index);
mask_manu = image_data.mask_manual;
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
fprintf('%s: pad image to accomodate edge displacements\n', mfilename);

pad_dim = ceil(args.pad_width*size(img));
img = padarray(img, pad_dim, 0, 'both');
roi = padarray(roi, pad_dim, 0, 'both');

% create new coordinate vectors, pixel units with origin at padded image center
[x_img, y_img] = meshgrid(1:size(img, 2), 1:size(img, 1));
x_img = x_img - mean(x_img(:));
y_img = y_img - mean(y_img(:));

% generate synthetic images  ----------

% compute exact displacement field for specified displacements and boundary
fprintf('%s: compute exact displacement field\n', mfilename);

[u_img, v_img] = exact_velocity(...
    x_img, y_img, args.defm_width/2*size(img, 2), args.theta, ...
    args.velocity_a, args.velocity_b);

fprintf('%s: generate synthetic images\n', mfilename);

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

% DEBUG: disabled plot
% % display initial and final synthetic images
% figure('Position', get(0, 'ScreenSize'))
% 
% subplot(1,2,1)
% imagesc([x_img(1), x_img(end)], [y_img(1), y_img(end)], ini);
% set(gca, 'YDir', 'normal', 'XGrid', 'on', 'YGrid', 'on', 'GridColor', 'w');
% axis equal tight
% title('Initial Synthetic Image')
% 
% subplot(1,2,2)
% imagesc([x_img(1), x_img(end)], [y_img(1), y_img(end)], fin);
% set(gca, 'YDir', 'normal', 'XGrid', 'on', 'YGrid', 'on', 'GridColor', 'w');
% axis equal tight
% title('Final Synthetic Image')

% run PIV analysis on synthetic images, get exact solution on same grid  ----------

fprintf('%s: run PIV analysis\n', mfilename);

result = piv_step(...     
    ini, fin, ini_roi, fin_roi, x_img(1,:), y_img(:,1), args.samp_len, args.samp_spc, ...
    args.intr_len, args.num_pass, args.valid_radius, args.valid_max, ...
    args.valid_eps, args.min_frac_data, args.min_frac_overlap);

% repackage key results
y_prime = result.x_grd*sind(args.theta) + result.y_grd*cosd(args.theta);
coord_obs = struct('x', result.x_grd, 'y', result.y_grd, 'roi', result.roi_grd, ...
    'y_prime', y_prime);
velocity_obs = struct('u', result.u_grd, 'v', result.v_grd, ...
    'm', sqrt(result.u_grd.^2 + result.v_grd.^2), ...
    'theta', atand(result.v_grd./result.u_grd));

% apply roi mask to velocity fields
velocity_obs.u(~coord_obs.roi) = NaN;
velocity_obs.v(~coord_obs.roi) = NaN;
velocity_obs.m(~coord_obs.roi) = NaN;

% compute exact velocities and package results
[u_ext, v_ext] = exact_velocity(...
    coord_obs.x, coord_obs.y, args.defm_width/2*size(img, 2), args.theta, ...
    args.velocity_a, args.velocity_b);
velocity_ext = struct('u', u_ext, 'v', v_ext, ...
    'm', sqrt(u_ext.^2 + v_ext.^2), ...
    'theta', atand(v_ext./u_ext));

% compute velocity errors and store in struct  ----------

velocity_error = struct();
velocity_fields = fieldnames(velocity_ext);
for ii = 1:numel(velocity_fields)
    fn = velocity_fields{ii};
    velocity_error.(fn) = velocity_ext.(fn) - velocity_obs.(fn);
end

% run deformation analysis, get exact deformation parameters on same grid  ----------

strain_obs = post_strain(coord_obs.x(1,:), coord_obs.y(:,1), velocity_obs.u, ...
                         velocity_obs.v, coord_obs.roi, 'nearest');

strain_ext = exact_strain(...
    coord_obs.x, coord_obs.y, args.defm_width/2*size(img, 2), args.theta, ...
    args.velocity_a, args.velocity_b);

% compute strain errors and store in struct ----------

strain_error = struct();
strain_fields = fieldnames(strain_ext);
for ii = 1:numel(strain_fields)
    fn = strain_fields{ii};
    strain_error.(fn) = strain_ext.(fn) - strain_obs.(fn);
end

% summarize error distributions ----------

fields = {'u', 'v', 'm', 'theta'};
names = {'u', 'v', 'velocity magnitude', '\theta'};
units = {'pixels/step', 'pixels/step', 'pixels/step', 'degrees'};
for ii = 1:length(names)
    plot_error_summary(coord_obs.x, coord_obs.y, coord_obs.y_prime, ...
        velocity_ext.(fields{ii}), velocity_obs.(fields{ii}), ...
        names{ii}, units{ii});
end

% TODO: add plots for deformation parameters

% save test inputs and outputs to file ----------

% TODO: make this an optional input

return

function plot_error_summary(xx, yy, zz, ext, obs, name, units)
% 
% Plot figure summarizing spatial distributions of absolute and relative errors
%
% Arguments: 
%   xx, yy: 2D coordinate matrices for error map
%   zz: 2D coordinate matrix for functional boxplot
%   ext: 2D matrix, Exact value of measured quantity
%   obs: 2D matrix, Observed value of measured quantity
%   name: String, variable name
%   units: String, variable units
% %

% define constants
screen_position = get(0, 'screensize');
fig_color = [1, 1, 1];
font_size = 10;

% create new figure
hf = figure;
hf.Color = fig_color;
hf.Position(1) = 1;
hf.Position(3) = screen_position(3); 

% compute errors
absolute = obs - ext;
relative = abs(absolute./ext);  % may contain inf

subplot(2, 2, 1);
plot_error_map(xx, yy, absolute, 'X Position [pixels]', 'Y Position [pixels]', ...
    sprintf('%s_{obs}-%s_{ext} [%s]', name, name, units), ...
    sprintf('%s Absolute Error Map', name));

subplot(2, 2, 2);
plot_error_dist(zz, absolute, 'Distance to Shear Zone Center [pixels]', ...
    sprintf('%s_{obs}-%s_{ext} [%s]', name, name, units), ...
    sprintf('%s Absolute Error Distribution', name));

subplot(2, 2, 3);
plot_error_map(xx, yy, relative, 'X Position [pixels]', 'Y Position [pixels]', ...
    sprintf('|%s_{obs}-%s_{ext}/%s_{ext}| [1]', name, name, name), ...
    sprintf('%s Relative Error Map', name));

subplot(2, 2, 4);
plot_error_dist(zz, relative, 'Distance to Shear Zone Center [pixels]', ...
    sprintf('|%s_{obs}-%s_{ext}/%s_{ext}| [1]', name, name, name), ...
    sprintf('%s Relative Error Distribution', name));

% TODO: style up the figure

return


function plot_error_map(xx, yy, ee, xlbl, ylbl, elbl, ttl)
% Create a spatial map of observed errors
%
% Arguments:
%   xx, yy: 2D coordinate matrices
%   ee: 2D error matrix, NaN where not observed
%   xlbl, ylbl: Strings, coordinate axis labels
%   elbl: String, colorbar label
%   ttl: String, plot title
% %
    imagesc([min(xx(:)), max(xx(:))], [min(yy(:)), max(yy(:))], ee, ...
        'AlphaData', ~isnan(ee));
    hcb = colorbar('EastOutside');
    xlabel(xlbl);
    ylabel(ylbl);
    hcb.Label.String = elbl;
    title(ttl);
    axis('equal', 'tight');
return


% TODO: cleanup botev code, make print statements optional
function plot_error_dist(zz, err, zlbl, elbl, ttl)
% Plot quantile ranges of the conditional distribution: P(err | zz)
%
% Arguments:
%   zz: 2D matrix, conditioning variable, typically some position value
%   err: 2D error matrix, NaN where not observed
%   zlbl, elbl: Strings, coordinate axis labels for zz and err respectively
%   ttl: String, plot title
% %

% define constants
e_num_grd = 1000;
z_num_grd = 100;
light_gray = 0.8;
dark_gray = 0.5;
fracpad = 0.3;
quant_rng = [0.99, 0.95, 0.5];

% crop data to roi
roi = ~isnan(err);
zz = zz(roi);
err = err(roi);

% make coords, padded edges allow for kernel smearing beyond data range
zpad = fracpad*range(zz(:)); 
zvec = linspace(min(zz(:)) - zpad, max(zz(:)) + zpad, z_num_grd);
epad = fracpad*range(err(:)); 
evec = linspace(min(err(:)) - epad, max(err(:)) + epad, e_num_grd);
espc = evec(2) - evec(1);
[zgrd, egrd] = meshgrid(zvec, evec);

% estimate distributions using adaptive kernel
joint_pdf = akde([zz(:), err(:)], [zgrd(:), egrd(:)], numel(zz));   
joint_pdf = reshape(joint_pdf, e_num_grd, z_num_grd);
marg_pdf = trapz(joint_pdf)*espc;
cond_pdf = bsxfun(@rdivide, joint_pdf, marg_pdf);
cond_cdf = cumtrapz(cond_pdf)*espc;

% interpolate quantiles
quant_bot = 0.5 - 0.5*quant_rng;
quant_top = 0.5 + 0.5*quant_rng;
cond_quant_bot = nan(numel(zvec), numel(quant_bot));
cond_quant_top = nan(numel(zvec), numel(quant_top));
cond_median = nan(numel(zvec), 1);
for ii = 1:numel(zvec)
    [uniq_cdf, uniq_idx] = unique(cond_cdf(:, ii)); % x must be unique to interp
    uniq_err = evec(uniq_idx);
    cond_quant_bot(ii, :) = interp1(uniq_cdf, uniq_err, quant_bot);
    cond_quant_top(ii, :) = interp1(uniq_cdf, uniq_err, quant_top);
    cond_median(ii) = interp1(uniq_cdf, uniq_err, 0.50);
end

% plot inter-quantile ranges
lgnd_txt = cell(numel(quant_rng) + 1 , 1);
grays = linspace(light_gray, dark_gray, 3)'*ones(1, 3); 
x_patch = [zvec, zvec(end:-1:1)];
for ii = 1:numel(quant_rng)
    lgnd_txt{ii} = sprintf('%.0f% (%.3f-%.3f)', 100*quant_rng(ii), ...
        quant_bot(ii), quant_top(ii));  
    y_patch = [cond_quant_top(:, ii); cond_quant_bot(end:-1:1, ii)]';
    patch(x_patch, y_patch, grays(ii,:), 'LineStyle', 'none');
end

% plot median
hold on
lgnd_txt{end} = 'Median';
plot(zvec, cond_median, 'k', 'LineWidth', 1); 

% format axes
legend(lgnd_txt);
xlabel(zlbl);
ylabel(elbl);
title(ttl);
set(gca, 'XLim', [min(zz(:)), max(zz(:))]);
grid on

return

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