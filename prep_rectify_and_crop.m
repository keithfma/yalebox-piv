%function [] = prep_rectify_and_crop(woco, img, )
%
% Use manually defined world-coordinate control points to project image into a
% regular coordinate system, and crop this image to the desired size.

load k24_woco_side.mat
im = imread('woco.jpg');
origin_xw = -0.6;
origin_yw = 0.0;
crop_xw = [-0.6, 0.7];
crop_yw = [0.0, 0.3];

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
ctrl_xr = 1+(ctrl_xw-origin_xw)*delta_pixel_per_meter;
ctrl_yr = 1+(ctrl_yw-origin_yw)*delta_pixel_per_meter;
 
% transform image to rectified coordinates
tform = fitgeotrans([ctrl_xp, ctrl_yp], [ctrl_xr, ctrl_yr],'lwm', 10);
warning('off', 'images:inv_lwm:cannotEvaluateTransfAtSomeOutputLocations');
imr = imwarp(im, tform);
warning('on', 'images:inv_lwm:cannotEvaluateTransfAtSomeOutputLocations');

% crop rectified image
crop_xr = 1+(crop_xw-origin_xw)*delta_pixel_per_meter;
crop_yr = 1+(crop_yw-origin_yw)*delta_pixel_per_meter;
cropbox = [crop_xr(1), crop_yr(1), diff(crop_xr), diff(crop_yr)];
imrc = imcrop(imr, cropbox);