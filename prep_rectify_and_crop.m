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

% option 1: get pixel size from pairwise distance between neighboring grid points
dist_world = [];
dist_pixel = [];
for ii = 1:num_ctrl
    
    nbrs = (abs(ctrl_xw(ii)-ctrl_xw) < (1.1*woco_world_dx)) & ...
        (abs(ctrl_yw(ii)-ctrl_yw) < (1.1*woco_world_dy)) & ...
        (ctrl_index ~= ii);
    
    nbrs_dist_world = sqrt((ctrl_xw(ii)-ctrl_xw(nbrs)).^2 + ...
        (ctrl_yw(ii)-ctrl_yw(nbrs)).^2);
    
    dist_world = [dist_world; nbrs_dist_world(:)]; %#ok!
    
    nbrs_dist_pixel = sqrt((ctrl_xp(ii)-ctrl_xp(nbrs)).^2 + ...
        (ctrl_yp(ii)-ctrl_yp(nbrs)).^2);
    
    dist_pixel = [dist_pixel; nbrs_dist_pixel(:)]; %#ok!
end

figure
hist(dist_world./dist_pixel*1000, 25);
title('neighboring points only');
fprintf('neighboring points only: delta = %.4f mm/pixel\n', median(dist_world./dist_pixel*1000));

% option 2: get pixel size from pairwise distance between all points
Dw = pdist([ctrl_xw, ctrl_yw]);
Dp = pdist([ctrl_xp, ctrl_yp]);

figure
hist(Dw./Dp*1000, 25);
title('all points');
fprintf('all points: delta = %.4f mm/pixel\n', median(Dw./Dp*1000));
