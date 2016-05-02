%function [] = prep_rectify_and_crop(woco)
%
% Use manually defined world-coordinate control points to project image into a
% regular coordinate system, and crop this image to the desired size.

load k24_woco_side.mat
woco_world_dx = 0.1; % used to find neighbors
woco_world_dy = 0.1;

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
