%function [rgb_r, coord_x, coord_y] = prep_rectify_and_crop(woco, rgb, crop_xlim_world, crop_ylim_world)
%
% Use manually defined world-coordinate control points to project image into a
% regular coordinate system, and crop this image to the desired size.

load k24_woco_side.mat
rgb = imread('woco.jpg');
crop_xlim_world = [-0.6, 0.6];
crop_ylim_world = [0.0, 0.25];

%  breakup woco struct
ctrl_xp = woco.xp;
ctrl_yp = woco.yp;
ctrl_xw = woco.xw;
ctrl_yw = woco.yw;
num_ctrl = numel(ctrl_xp);
ctrl_index = (1:num_ctrl)';

% get pixel size from median of pairwise distance between all points
Dw = pdist([ctrl_xw, ctrl_yw]);
Dp = pdist([ctrl_xp, ctrl_yp]);
delta_pixel_per_meter = median(Dp./Dw);
delta_meter_per_pixel = median(Dw./Dp);

% get rectified pixel coordinates from origin and pixel size
ctrl_xr = 1+(ctrl_xw-crop_xlim_world(1))*delta_pixel_per_meter;
ctrl_yr = 1+(ctrl_yw-crop_ylim_world(1))*delta_pixel_per_meter;
 
% transform image to rectified coordinates
warning('off', 'images:inv_lwm:cannotEvaluateTransfAtSomeOutputLocations');
tform = fitgeotrans([ctrl_xp, ctrl_yp], [ctrl_xr, ctrl_yr],'lwm', 10);
rgb_r = imwarp(rgb, tform);
warning('on', 'images:inv_lwm:cannotEvaluateTransfAtSomeOutputLocations');

% crop rectified image
crop_xlim_r = 1+(crop_xlim_world-crop_xlim_world(1))*delta_pixel_per_meter;
crop_ylim_r = 1+(crop_ylim_world-crop_ylim_world(1))*delta_pixel_per_meter;
cropbox = [crop_xlim_r(1), crop_ylim_r(1), diff(crop_xlim_r), diff(crop_ylim_r)];
rgb_r = imcrop(rgb_r, cropbox);

% generate coordinate vectors (assumes that x=0, y = 0 is a control point)
origin_xr = ctrl_xr(ctrl_xw == 0  & ctrl_yw == 0);
origin_yr = ctrl_yr(ctrl_xw == 0  & ctrl_yw == 0);
coord_x = ((1:size(rgb_r, 2))-origin_xr)*delta_meter_per_pixel;
coord_y = ((1:size(rgb_r, 1))-origin_yr)*delta_meter_per_pixel;