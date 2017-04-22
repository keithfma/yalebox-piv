function [] = test_piv(varargin)
% function [] = test_piv(varargin)
%
% Test PIV performance using case using synthetic image pair generated by
% homogenous deformation and translation of a single input image.
%
% NOTE: This is a "high-level" test that covers all of the PIV subsystems. The
% error analysis addresses only measured velocities, not derived products such
% as strain.
%
% Optional Arguments ('Name', Value):
%   'image_file': String, name of pre-processed image file, as produced by
%       prep_series(), to read for raw image, default = ./test/default_image.nc
%   'image_index': Integer, (1-based) index of image in image file to use for
%       raw image, default = 1
%   'image_pos': 4-element position vector indicating the limits of the image to
%       extract and deform. Must contain only sand (all within the ROI),
%       in meters, default = [-0.12, 0.005, 0.092, 0.07]
%   'translation': 2-element vector specifying spatially constant translation
%       in meters, default = [0.005, 0.00]
%   'shear_theta': Scalar, orientation of shear band specified as
%       counter-clockwise angle to the positive x-axis, in degrees, limited to
%       range 0 - 90, default = 45
%   'shear_width': Scalar, width of shear band, in meters, default = 0.05   
%   'shear_mag': Scalar, displacement difference across shear band, applied as a
%       0.5*shear_mag displacement on one side and a -0.5*shear_mag displacement
%       on the other, default = sqrt(2)*0.005
%   'pad_width': Scalar, width of edge padding to add to image (to accomodate
%       edge displacements) as a fraction of image size, default = 0.1
%   'samplen': piv() parameter, default [30, 30]
%   'sampspc': piv() parameter, default [15, 15]
%   'intrlen': piv() parameter, default [100, 60]
%   'npass': piv() parameter, default [1, 2]
%   'valid_max': piv() parameter, default 2
%   'valid_eps': piv() parameter, default 0.1
%   'spline_tension': piv() parameter, default 0.95
%   'min_frac_data': piv() parameter, default 0.8
%   'min_frac_overlap': piv() parameter, default 0.5
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
ip.addParameter('translation', [0.005, 0.00], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('shear_theta', 45, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=' 90}));
ip.addParameter('shear_width', 0.05, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
ip.addParameter('shear_mag', 0.01, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
ip.addParameter('pad_width', 0.1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0}));
ip.addParameter('samplen', [30, 30]); % validation handled by PIV routines
ip.addParameter('sampspc', [15, 15]);
ip.addParameter('intrlen', [100, 60]);
ip.addParameter('npass', [1, 2]);
ip.addParameter('valid_max', 2);
ip.addParameter('valid_eps', 0.1);
ip.addParameter('spline_tension', 0.95);
ip.addParameter('min_frac_data', 0.8);
ip.addParameter('min_frac_overlap', 0.5);
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

xx = double(ncread(args.image_file, 'x'));
yy = double(ncread(args.image_file, 'y'));
img = double(ncread(args.image_file, 'img', [1, 1, args.image_index], [inf, inf, 1]));
mask_auto = ncread(args.image_file, 'mask_auto', [1, 1, args.image_index], [inf, inf, 1]);
mask_manu = ncread(args.image_file, 'mask_manual');
roi = mask_auto & mask_manu;

min_col = find(xx >= args.image_pos(1), 1, 'first');
max_col = find(xx <= args.image_pos(1) + args.image_pos(3), 1, 'last');
min_row = find(yy >= args.image_pos(2), 1, 'first');
max_row = find(yy <= args.image_pos(2) + args.image_pos(4), 1, 'last');

xx = xx(min_col:max_col); % not used, included for reference
yy = yy(min_row:max_row);
img = img(min_row:max_row, min_col:max_col);
roi = roi(min_row:max_row, min_col:max_col);

if any(~roi(:))
    error('%s: image_pos limits must include only sand (ROI)', mfilename);
end

%% convert units from world (meters) to image (pixels)
% NOTE: all the analysis from here on out is in pixels, these are the units that
%   matter to the PIV routines.

% conversion factors
m_per_pix = mean([diff(xx); diff(yy)]); % grid is regular to single precision
pix_per_m = 1/m_per_pix;

% convert parameters
args.translation = args.translation*pix_per_m;
args.shear_width = args.shear_width*pix_per_m;
args.shear_mag = args.shear_mag*pix_per_m;

%% pad image boundaries (and coordinates) to accomodate edge displacements

if args.verbose
    fprintf('%s: pad image to accomodate edge displacements\n', mfilename);
end

pad_dim = ceil(args.pad_width*size(img));
img = padarray(img, pad_dim, 0, 'both');
roi = padarray(roi, pad_dim, 0, 'both');

% create new coordinate vectors, pixel units with origin at padded image center
x_img = (1:size(img, 2))-1;
y_img = (1:size(img, 1))-1;
x_img = x_img - mean(x_img);
y_img = y_img - mean(y_img);
[x_img_grd, y_img_grd] = meshgrid(x_img, y_img);

%% compute exact displacement field for specified displacements and boundary

if args.verbose
    fprintf('%s: compute exact displacement field\n', mfilename);
end

[u_exact_img, v_exact_img] = get_exact_uv(...
    x_img_grd, y_img_grd, args.translation, args.shear_theta, ...
    args.shear_width, args.shear_mag);

%% generate synthetic images

if args.verbose
    fprintf('%s: generate synthetic images\n', mfilename);
end

% deform images
ini = imwarp(img, 0.5*cat(3, u_exact_img, v_exact_img), 'cubic'); % dir is ok
fin = imwarp(img, -0.5*cat(3, u_exact_img, v_exact_img), 'cubic');

% deform masks
ini_roi = imwarp(double(roi), 0.5*cat(3, u_exact_img, v_exact_img), 'cubic');
ini_roi = logical(round(ini_roi));
fin_roi = imwarp(double(roi), -0.5*cat(3, u_exact_img, v_exact_img), 'cubic');
fin_roi = logical(round(fin_roi));

% reapply mask
ini(~ini_roi) = 0;
fin(~fin_roi) = 0;

% enforce limits by threshold
% NOTE: could stretch limits instead, unclear if this matters
ini(ini < 0) = 0;
ini(ini > 1) = 1;
fin(fin < 0) = 0;
fin(fin > 1) = 1;

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

%% run PIV analysis on synthetic images

if args.verbose
    fprintf('%s: run PIV analysis\n', mfilename);
end

[x_piv, y_piv, u_piv, v_piv, roi_piv] = piv(... 
    ini, fin, ini_roi, fin_roi, x_img, y_img, args.samplen, args.sampspc, ...
    args.intrlen, args.npass, args.valid_max, args.valid_eps, ...
    args.spline_tension, args.min_frac_data, args.min_frac_overlap, ...
    args.verbose);

% display exact and measure displacement fields
figure('Position', get(0, 'ScreenSize'))

num_vec = 25; % desired num quiver vectors along largest dim, for downsampling

subplot(1, 2, 1)
dfact = floor(max(size(u_exact_img))/num_vec);
m_exact = sqrt(u_exact_img.^2 + v_exact_img.^2);
imagesc(x_img, y_img, m_exact);
set(gca, 'YDir', 'Normal', 'XGrid', 'on', 'YGrid', 'on', 'GridColor', 'w');
hold on;
quiver(x_img_grd(1:dfact:end, 1:dfact:end), ...
    y_img_grd(1:dfact:end, 1:dfact:end), ...
    u_exact_img(1:dfact:end, 1:dfact:end), ...
    v_exact_img(1:dfact:end, 1:dfact:end), '-k');
cb = colorbar;
cb.Label.String = 'Displacement Magnitude [pixels]';
axis equal tight
title('Exact Displacement Magnitude and Direction')

subplot(1, 2, 2)
dfact = floor(max(size(u_piv))/num_vec);
m_piv = sqrt(u_piv.^2 + v_piv.^2);
imagesc(x_piv, y_piv, m_piv, 'AlphaData', roi_piv);
set(gca, 'YDir', 'Normal', 'XGrid', 'on', 'YGrid', 'on', 'GridColor', 'w');
hold on;
[x_piv_grd, y_piv_grd] = meshgrid(x_piv, y_piv);
quiver(x_piv_grd(1:dfact:end, 1:dfact:end), ...
    y_piv_grd(1:dfact:end, 1:dfact:end), ...
    u_piv(1:dfact:end, 1:dfact:end), ...
    v_piv(1:dfact:end, 1:dfact:end), '-k');
cb = colorbar;
cb.Label.String = 'Displacement Magnitude [pixels]';
axis equal tight
title('PIV Displacement Magnitude and Direction')


%% analyze errors

if args.verbose
    fprintf('%s: error analysis\n', mfilename);
end

[u_exact_at_piv, v_exact_at_piv] = get_exact_uv(...
    x_piv_grd, y_piv_grd, args.translation, args.shear_theta, ...
    args.shear_width, args.shear_mag);
u_error = u_exact_at_piv - u_piv;
v_error = v_exact_at_piv - v_piv;
m_error = sqrt(u_error.^2 + v_error.^2);
theta_error = atan2d(v_exact_at_piv, u_exact_at_piv) - atan2d(v_piv, u_piv);

% print error quantiles
qnt = 0 : 0.10 : 1;
u_abs_error_qnt = quantile(abs(u_error(:)), qnt);
v_abs_error_qnt = quantile(abs(v_error(:)), qnt);
m_error_qnt = quantile(sqrt(u_error(:).^2 + v_error(:).^2), qnt);

fprintf('\n---- \n');
fprintf('Displacement Vector Error Quantiles\n\n');
fprintf('Quantile\t| U Abs Error\t| V Abs Error\t| Mag Error\n'); 
for i = 1:length(qnt)
    fprintf('%.2f\t| %.2e\t| %.2e\t| %.2e\n', ...
        qnt(i), u_abs_error_qnt(i), v_abs_error_qnt(i), m_error_qnt(i));
end
fprintf('---- \n\n');

% plot error maps
figure('Position', get(0, 'ScreenSize'));

subplot(2, 2, 1);
imagesc(x_piv, y_piv, u_error, 'AlphaData', roi_piv);
set(gca, 'YDir', 'normal');
axis equal tight
cb = colorbar;
cb.Label.String = 'Error [pixels]';
title('U Error')

subplot(2, 2, 2);
imagesc(x_piv, y_piv, v_error, 'AlphaData', roi_piv);
set(gca, 'YDir', 'normal');
axis equal tight
cb = colorbar;
cb.Label.String = 'Error [pixels]';
title('V Error')

subplot(2, 2, 3);
imagesc(x_piv, y_piv, m_error, 'AlphaData', roi_piv);
set(gca, 'YDir', 'normal');
axis equal tight
cb = colorbar;
cb.Label.String = 'Error [pixels]';
title('Magnitude Error')

subplot(2, 2, 4);
imagesc(x_piv, y_piv, theta_error, 'AlphaData', roi_piv);
set(gca, 'YDir', 'normal');
axis equal tight
cb = colorbar;
cb.Label.String = 'Error [degrees]';
title('Theta Error')

% plot error histograms
num_bins = 50;

figure('Position', get(0, 'ScreenSize'))

subplot(2, 2, 1)
hist(u_error(~isnan(u_error)), num_bins);
set(gca, 'XLim', [min(u_error(:)), max(u_error(:))]);
hh = findobj(gca, 'Type', 'patch');
hh.LineStyle = 'none';
title('U Error')
xlabel('Error [pixels]');
ylabel('Count [pixels]');

subplot(2, 2, 2)
hist(v_error(~isnan(v_error)), num_bins);
set(gca, 'XLim', [min(v_error(:)), max(v_error(:))]);
hh = findobj(gca, 'Type', 'patch');
hh.LineStyle = 'none';
title('V Error')
xlabel('Error [pixels]');
ylabel('Count [pixels]');

subplot(2, 2, 3)
hist(m_error(~isnan(m_error)), num_bins);
set(gca, 'XLim', [min(m_error(:)), max(m_error(:))]);
hh = findobj(gca, 'Type', 'patch');
hh.LineStyle = 'none';
title('Magnitude Error')
xlabel('Error [pixels]');
ylabel('Count [pixels]');

subplot(2, 2, 4)
hist(theta_error(~isnan(theta_error)), num_bins);
set(gca, 'XLim', [min(theta_error(:)), max(theta_error(:))]);
hh = findobj(gca, 'Type', 'patch');
hh.LineStyle = 'none';
title('Theta Error')
xlabel('Error [degrees]');
ylabel('Count [pixels]');


function [u_grd, v_grd] = get_exact_uv(...
        x_grd, y_grd, translation, shear_theta, shear_width, shear_mag)
% function [uu, vv] = get_exact_uv(...
%         xx, yy, translation, shear_theta, shear_width, shear_mag)
%
% Compute exact displacement field at specified coordinate grid
%
% Arguments:
%   x_grd, y_grd: Coordinate grids, as produced by meshgrid
%   all others: displacement field parameters, see test_piv() help
% %

sz = size(x_grd);
u_grd = zeros(sz);
v_grd = zeros(sz);

% apply constant displacement
u_grd = u_grd + translation(1);
v_grd = v_grd + translation(2);

% create rotated coordinate system with y == 0 at center of shear band
rot = [cosd(shear_theta), sind(shear_theta); ...
       -sind(shear_theta), cosd(shear_theta)];
xy_rot = rot*[x_grd(:)'; y_grd(:)'];
y_grd_rot = reshape(xy_rot(2,:), sz);

% compute scaling factor for shear displacements
shear_scale = y_grd_rot/shear_width;
shear_scale(shear_scale < -0.5) = -0.5;
shear_scale(shear_scale > 0.5) = 0.5;

% apply shear displacements
u_shear = shear_scale*cosd(shear_theta)*shear_mag;
v_shear = shear_scale*sind(shear_theta)*shear_mag;
u_grd = u_grd + u_shear;
v_grd = v_grd + v_shear;
