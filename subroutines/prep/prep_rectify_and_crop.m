function [rgb_out, coord_x, coord_y] = prep_rectify_and_crop(...
    ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, rgb, show)
% function [rgb_out, coord_x, coord_y] = prep_rectify_and_crop(...
%     ctrl_xp, ctrl_yp, ctrl_xw, ctrl_yw, crop_xw, crop_yw, rgb, show)
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
% rgb_r = 3D matrix, rectified and cropped image
%
% coord_x, coord_y = Vectors, world coordinate vectors for rgb_r
%
% %
warning('off', 'images:inv_lwm:cannotEvaluateTransfAtSomeOutputLocations');
warning('off', 'images:geotrans:estimateOutputBoundsFailed');

% constants
FIT_NPTS = 10; 

% set defaults
if nargin < 8 || isempty(show); show = false; end

% sanity check
validateattributes(ctrl_xp, {'numeric'}, {'vector'});
validateattributes(ctrl_yp, {'numeric'}, {'vector', 'numel', numel(ctrl_xp)});
validateattributes(ctrl_xw, {'numeric'}, {'vector', 'numel', numel(ctrl_xp)});
validateattributes(ctrl_yw, {'numeric'}, {'vector', 'numel', numel(ctrl_xp)});
validateattributes(crop_xw, {'numeric'}, {'vector', 'numel', 2});
validateattributes(crop_yw, {'numeric'}, {'vector', 'numel', 2});
assert(crop_xw(2)>crop_xw(1));
assert(crop_yw(2)>crop_yw(1));

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
pixels_per_meter = median(Dp./Dw);
meters_per_pixel = 1/pixels_per_meter;

% adjust crop limits to register pixels to origin, i.e., ensure we have a pixel at
%   world coordinates (0, 0)
crop_xp = [floor(crop_xw(1)*pixels_per_meter) - 0.5, ceil(crop_xw(2)*pixels_per_meter) + 0.5];
crop_yp = [floor(crop_yw(1)*pixels_per_meter) - 0.5, ceil(crop_yw(2)*pixels_per_meter) + 0.5];
crop_xw = crop_xp * meters_per_pixel;
crop_yw = crop_yp * meters_per_pixel;

% create output image reference
out_size = [diff(crop_yp), diff(crop_xp)];   % plus one?
out_ref = imref2d(out_size, crop_xw, crop_yw);

% transform image to rectified pixel coordinates
tform = fitgeotrans([ctrl_xp, ctrl_yp], [ctrl_xw, ctrl_yw], 'lwm', FIT_NPTS);
in_ref = imref2d([size(rgb,1), size(rgb,2)]);
rgb_out = imwarp(rgb, in_ref, tform, 'cubic', 'OutputView', out_ref);

% generate world coordinate vectors for rectified, cropped image
[coord_x0, coord_y0] = out_ref.intrinsicToWorld(0, 0);  % intentionally not 1,1
coord_x = coord_x0 + out_ref.PixelExtentInWorldX*(1:size(rgb_out, 2));
coord_y = coord_y0 + out_ref.PixelExtentInWorldY*(1:size(rgb_out, 1));

fprintf('%s: output image size = %d x %d x %d\n', ...
    mfilename, size(rgb_out, 1), size(rgb_out, 2), size(rgb_out, 3));

% optional: show results
if show
   hf = figure;
   hf.Name = 'prep_rectify_and_crop: results';
   imagesc(coord_x, coord_y, rgb_out)
   set(gca, 'YDir', 'normal');
   daspect([1, 1, 1]);
   xlabel('x-position [m]');
   ylabel('y-position [m]');
end