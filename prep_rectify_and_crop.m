function [rgb_r, coord_x, coord_y] = ...
    prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, rgb, show, verbose)
% function [rgb_r, coord_x, coord_y] = ...
%     prep_rectify_and_crop(ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, rgb, show, verbose)
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
% show = Logical flag, plot the rectified image, default = false 
%
% show = Logical flag, print verbose messages, default = false 
%
% rgb_r = 3D matrix, rectified and cropped image
%
% coord_x, coord_y = Vectors, world coordinate vectors for rgb_r
%
% % Keith Ma

% sanity check
if nargin == 7 || isempty(show); show = false; end
if nargin == 8 || isempty(verbose); verbose = false; end

validateattributes(ctrl_xp, {'numeric'}, {'vector'});
validateattributes(ctrl_yp, {'numeric'}, {'vector', 'numel', numel(ctrl_xp)});
validateattributes(ctrl_xw, {'numeric'}, {'vector', 'numel', numel(ctrl_xp)});
validateattributes(ctrl_yw, {'numeric'}, {'vector', 'numel', numel(ctrl_xp)});
validateattributes(crop_xw, {'numeric'}, {'vector', 'numel', 2});
validateattributes(crop_yw, {'numeric'}, {'vector', 'numel', 2});
assert(crop_xw(2)>crop_xw(1));
assert(crop_yw(2)>crop_yw(1));

if verbose
    fprintf('%s: input image size = %d x %d x %d\n', ...
        mfilename, size(rgb, 1), size(rgb, 2), size(rgb, 3));
end

% get pixel size from median of pairwise distance between all points
Dw = pdist([ctrl_xw, ctrl_yw]);
Dp = pdist([ctrl_xp, ctrl_yp]);
delta_pixel_per_meter = median(Dp./Dw);
delta_meter_per_pixel = median(Dw./Dp);

% get rectified pixel coordinates from origin and pixel size
ctrl_xr = 1+(ctrl_xw-crop_xw(1))*delta_pixel_per_meter;
ctrl_yr = 1+(ctrl_yw-crop_yw(1))*delta_pixel_per_meter;
 
% transform image to rectified coordinates
warning('off', 'images:inv_lwm:cannotEvaluateTransfAtSomeOutputLocations');
warning('off', 'images:geotrans:estimateOutputBoundsFailed');
tform = fitgeotrans([ctrl_xp, ctrl_yp], [ctrl_xr, ctrl_yr],'lwm', 10);
rgb_r = imwarp(rgb, tform);
warning('on', 'images:inv_lwm:cannotEvaluateTransfAtSomeOutputLocations');
warning('on', 'images:geotrans:estimateOutputBoundsFailed');

% crop rectified image
crop_xlim_r = 1+(crop_xw-crop_xw(1))*delta_pixel_per_meter;
crop_ylim_r = 1+(crop_yw-crop_yw(1))*delta_pixel_per_meter;
cropbox = [crop_xlim_r(1), crop_ylim_r(1), diff(crop_xlim_r), diff(crop_ylim_r)];
rgb_r = imcrop(rgb_r, cropbox);

% generate coordinate vectors (assumes that x=0, y = 0 is a control point)
origin_xr = ctrl_xr(ctrl_xw == 0  & ctrl_yw == 0);
origin_yr = ctrl_yr(ctrl_xw == 0  & ctrl_yw == 0);
coord_x = ((1:size(rgb_r, 2))-origin_xr)*delta_meter_per_pixel;
coord_y = ((1:size(rgb_r, 1))-origin_yr)*delta_meter_per_pixel;

if verbose
    fprintf('%s: output image size = %d x %d x %d\n', ...
        mfilename, size(rgb_r, 1), size(rgb_r, 2), size(rgb_r, 3));
end

% optional: show results
if show
   hf = figure;
   hf.Name = 'prep_rectify_and_crop: results';
   imagesc(coord_x, coord_y, rgb_r)
   set(gca, 'YDir', 'normal');
   xlabel('x-position [m]');
   ylabel('y-position [m]');
end