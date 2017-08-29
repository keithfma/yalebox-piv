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
%       default = 20*[cosd(45), sind(45)]
%   'u2': 2-element vector specifying second end-member translation in pixels,
%       default = 10*[cosd(45), sind(45)]
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
ip.addParameter('velocity_a', 10*[cosd(45), sind(45)], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('velocity_b', 20*[cosd(45), sind(45)], ...
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
[xx, yy] = meshgrid(1:size(img, 2), 1:size(img, 1));
xx = xx - mean(xx(:));
yy = yy - mean(yy(:));

if any(~roi(:))
    error('%s: image_pos limits must include only sand (ROI)', mfilename);
end

%% DEBUG

[uu, vv] = exact_velocity(...
    xx, yy, args.defm_width/2*size(img, 2), args.theta, args.velocity_a, ...
    args.velocity_b);

[dudx, dudy, dvdx, dvdy] = exact_gradient(...
    xx, yy, args.defm_width/2*size(img, 2), args.theta, args.velocity_a, ...
    args.velocity_b);

% TODO: compute exact strain parameters

keyboard

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

dudx(defm) = (uv_b(1) - uv_a(1))*sind(theta)/(2*hh);
dvdx(defm) = (uv_b(2) - uv_a(2))*sind(theta)/(2*hh);
dudy(defm) = (uv_b(1) - uv_a(1))*cosd(theta)/(2*hh);
dvdy(defm) = (uv_b(2) - uv_a(2))*cosd(theta)/(2*hh);

return

% TODO: compute strain parameters
