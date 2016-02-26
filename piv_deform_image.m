function [img_tm] = piv_deform_image(img_tx, img_roi_tx, r_grd_tm, c_grd_tm, u_grd_tm, v_grd_tm, roi, is_fwd)
%
% img = Original image
% r_pts, c_pts = Scattered grid for displacement estimates at image time
% u_pts, v_pts = Displacements on the scattered grid at image time
% is_fwd = Logical flag, deform the image a half-step forward in time (1), or a half step back in time (0).
% %

% tx = initial or final time, depending

% constant parameters
spc = 5; % low-res grid spacing in pixels
tension = 0.95;
roi_epsilon = 1e-2; % numerical precision for roi deformation

% get full-res grid of images coordinates
[nr_img, nc_img] = size(img_tx);
[c_img_tx, r_img_tx] = meshgrid(1:nc_img, 1:nr_img); % full-res

% get low-res grid that spans the image coordinates
cc = 1:spc:(ceil(nc_img/spc)*spc+1);
rr = 1:spc:(ceil(nr_img/spc)*spc+1);
[c_lr_tx, r_lr_tx] = meshgrid(cc, rr); 
[nr_lr, nc_lr] = size(c_lr_tx);

% propagate points to target time (half-step forward or back, depending)
if is_fwd
    % forward image defm, propagate displacements back from midpoint to initial time
    c_pts_tx = c_grd_tm - 0.5*u_grd_tm;
    r_pts_tx = r_grd_tm - 0.5*v_grd_tm;    
else
    % backward image defm, propagate displacements fwd from midpoint to final time
    c_pts_tx = c_grd_tm + 0.5*u_grd_tm;
    r_pts_tx = r_grd_tm + 0.5*v_grd_tm;    
end

% interpolate scattered to low-res grid using expensive tension splines
u_lr_tx = spline2d(c_lr_tx(:), r_lr_tx(:), c_pts_tx(roi), r_pts_tx(roi), ...
    u_grd_tm(roi), tension);
v_lr_tx = spline2d(c_lr_tx(:), r_lr_tx(:), c_pts_tx(roi), r_pts_tx(roi), ...
    v_grd_tm(roi), tension);
u_lr_tx = reshape(u_lr_tx, nr_lr, nc_lr);
v_lr_tx = reshape(v_lr_tx, nr_lr, nc_lr);

% interpolate on low-res to full-res grid using cheap linear interpolation
interp = griddedInterpolant(r_lr_tx, c_lr_tx, u_lr_tx, 'linear');
u_img_tx = interp(r_img_tx, c_img_tx);
interp.Values = v_lr_tx;
v_img_tx = interp(r_img_tx, c_img_tx);

% prep displacement matrix for image deformation 
if is_fwd
    displacement = -0.5*cat(3, u_img_tx, v_img_tx);
else
    displacement = 0.5*cat(3, u_img_tx, v_img_tx);
end

% deform image to midpoint time
img_tm = imwarp(img_tx, displacement, 'cubic', 'FillValues', 0);

% deform roi to midpoint time
tmp = imwarp(double(img_roi_tx), displacement, 'cubic', 'FillValues', 0);
img_roi_tm = abs(tmp-1) < roi_epsilon;

% re-apply roi
img_tm(~img_roi_tm) = 0;