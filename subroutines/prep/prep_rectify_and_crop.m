function [rgb_r, coord_x, coord_y] = ...
    prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, ...
                          rgb, fit_npts, show)
% function [rgb_r, coord_x, coord_y] = ...
%     prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, ...
%                           rgb, fit_npts, show)
% 
% Use manually defined world-coordinate control points to project image
% into a regular coordinate system, and crop this image to the desired
% size. The transformation is performed using the MATLAB image processing
% toolbox tools fitgeotrans and imwarp, and uses a local weighted mean
% transformation with cubic interpolation. This is a very flexible
% transform, and it is able to deal with lens distortion and perspective.
% 
% Arguments:
% 
% ctrl_xp, ctrl_yp = Vectors, control point location in pixel coordinates,
%   these are points in the rgb image
%
% ctrl_xw, ctrl_yw = Vectors, control point location in world coordinates,
%   these are points on the coordinate reference grid
%
% crop_xw, crop_yw = Vectors, numel == 2, [min, max] world coordinates for
%   the cropped output image
%
% rgb = 3D matrix, rgb image to be rectified and cropped
%
% fit_npts = Integer, parameter to fitgeotrans "Number of points to use in local
%   weighted mean calculation". The default value usually works well, but
%   occasionally leads to ill-conditioned arrays (throws a misleading "co-linear
%   points" error), default = 10
%
% show = Logical flag, plot the rectified image, default = false 
%
% rgb_r = 3D matrix, rectified and cropped image
%
% coord_x, coord_y = Vectors, world coordinate vectors for rgb_r
%
% %

% sanity check
if nargin < 8 || isempty(fit_npts); fit_npts = 10; end
if nargin < 9 || isempty(show); show = false; end

validateattributes(ctrl_xp, {'numeric'}, {'vector'});
validateattributes(ctrl_yp, {'numeric'}, {'vector', 'numel', numel(ctrl_xp)});
validateattributes(ctrl_xw, {'numeric'}, {'vector', 'numel', numel(ctrl_xp)});
validateattributes(ctrl_yw, {'numeric'}, {'vector', 'numel', numel(ctrl_xp)});
validateattributes(crop_xw, {'numeric'}, {'vector', 'numel', 2});
validateattributes(crop_yw, {'numeric'}, {'vector', 'numel', 2});
assert(crop_xw(2)>crop_xw(1));
assert(crop_yw(2)>crop_yw(1));
validateattributes(fit_npts', {'numeric'}, {'scalar', 'integer', 'positive'});

fprintf('%s: input image size = %d x %d x %d\n', ...
    mfilename, size(rgb, 1), size(rgb, 2), size(rgb, 3));

% ensure control points are column vectors
ctrl_xw = ctrl_xw(:);
ctrl_yw = ctrl_yw(:);
ctrl_xp = ctrl_xp(:);
ctrl_yp = ctrl_yp(:);

% get pixel size from median of pairwise distance between all points
Dw = pdist([ctrl_xw, ctrl_yw]);
Dp = pdist([ctrl_xp, ctrl_yp]);
delta_pixel_per_meter = median(Dp./Dw);
delta_meter_per_pixel = median(Dw./Dp);

% convert parameters to rectified pixel coordinates 
% ...origin at world (0,0), regular pixel size
ctrl_xr = ctrl_xw*delta_pixel_per_meter;
ctrl_yr = ctrl_yw*delta_pixel_per_meter;
crop_xr = crop_xw*delta_pixel_per_meter;
crop_yr = crop_yw*delta_pixel_per_meter;
origin_xr = 0;
origin_yr = 0;

% transform image to rectified pixel coordinates
warning('off', 'images:inv_lwm:cannotEvaluateTransfAtSomeOutputLocations');
warning('off', 'images:geotrans:estimateOutputBoundsFailed');
tform = fitgeotrans([ctrl_xp, ctrl_yp], [ctrl_xr, ctrl_yr],'lwm', fit_npts);
imref = imref2d([size(rgb,1), size(rgb,2)]);
[rgb_r, imref_r] = imwarp(rgb, imref, tform, 'cubic');
warning('on', 'images:inv_lwm:cannotEvaluateTransfAtSomeOutputLocations');
warning('on', 'images:geotrans:estimateOutputBoundsFailed');

% get subscripts for origin and crop limits in rectified image
[origin_col, origin_row] = imref_r.worldToIntrinsic(origin_xr, origin_yr);

[crop_col, crop_row] = imref_r.worldToIntrinsic(crop_xr, crop_yr);

crop_row(1) = max(1, floor(crop_row(1)));
crop_row(2) = min(size(rgb_r,1), ceil(crop_row(2)));

crop_col(1) = max(1, floor(crop_col(1)));
crop_col(2) = min(size(rgb_r,2), ceil(crop_col(2)));

crop_rect = [crop_col(1), crop_row(1), diff(crop_col), diff(crop_row)];

% crop rectified image, adjust origin location
rgb_r = imcrop(rgb_r, crop_rect);
origin_row = 1+origin_row-crop_row(1);
origin_col = 1+origin_col-crop_col(1);

% generate world coordinate vectors for rectified, cropped image
% ... convert pixel distance to origin to meters
coord_x = ((1:size(rgb_r, 2))-origin_col)*delta_meter_per_pixel;
coord_y = ((1:size(rgb_r, 1))-origin_row)*delta_meter_per_pixel;

fprintf('%s: output image size = %d x %d x %d\n', ...
    mfilename, size(rgb_r, 1), size(rgb_r, 2), size(rgb_r, 3));

% optional: show results
if show
   hf = figure;
   hf.Name = 'prep_rectify_and_crop: results';
   imagesc(coord_x, coord_y, rgb_r)
   set(gca, 'YDir', 'normal');
   daspect([1, 1, 1]);
   xlabel('x-position [m]');
   ylabel('y-position [m]');
end