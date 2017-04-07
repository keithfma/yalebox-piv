function img_tm = piv_deform_image(img_tx, img_roi_tx, r_grd_tm, c_grd_tm, ...
                      u_grd_tm, v_grd_tm, roi, tension, is_fwd, verbose)
% function img_tm = piv_deform_image(img_tx, img_roi_tx, r_grd_tm, c_grd_tm, ...
%                       u_grd_tm, v_grd_tm, roi, tension, is_fwd, verbose)
%
% Deform initial or final image to midpoint time based on current estimates for
% displacement. Initial displacement estimates are on a regular grid at midpoint
% time are propagated to the image time, then regridded at image resolution.
% These derived displacement estimates, at image resolution and image time, are
% used to deform the image forward or back one half time step.
%
% Arguments:
%
%   img_tx = 2D matrix, original image
%
%   img_roi_tx = 2D matrix, logical, indicates whether pixels have sand (1) or 
%       not (0)
%
%   r_grd_tm, c_grd_tm = 2D matrix, sample coordinate grid for displacement 
%       estimates at midpoint time
%
%   u_grd_tm, v_grd_tm = 2D matrix, estimated displacements on the regular 
%       sample grid at midpoint time
%
%   roi = 2D matrix, logical, indicates whether the PIV had enough data to 
%       estimate displacements for this point (1) or not (0)
%
%   tension = Scalar, tension parameter for the spline interpolation routine
%
%   is_fwd = Logical flag, deform the image a half-step forward in time (1), 
%       or a half step back in time (0).
%
%   verbose = Logical flag, display verbose messages (1) or don't
% %

% TODO: think about what interpolation method to use here...

% Note: A few meaningful suffixes are used to help clarify the variable grids
% and times, specifically:
% 
%   tm -> midpoint time
%   tx -> image time, which is either initial or final time depending on is_fwd
%   grd -> regular sample coordinate grid
%   img -> regular image coordinate grid
%   pts -> scattered coordinate points
%   lr -> regular "low-res" coordinate grid used in interpolation

% constants
roi_epsilon = 1e-2; % numerical threshold for roi deformation

if verbose
    fprintf('%s: tension = %.2f\n', mfilename, tension);
end

% get full-res grid of images coordinates
[nr_img, nc_img] = size(img_tx);
[c_img, r_img] = meshgrid(1:nc_img, 1:nr_img); % full-res

% get extended grid that spans the image coordinates
r_spc = r_grd_tm(2,1)-r_grd_tm(1,1);
r0 = r_grd_tm(1,1);
while r0 > 1
    r0 = r0-r_spc;
end
r1 = r_grd_tm(end,1);
while r1 < nr_img
    r1 = r1+r_spc;
end

c_spc = c_grd_tm(1,2)-c_grd_tm(1,1);
c0 = c_grd_tm(1,1);
while c0 > 1
    c0 = c0-c_spc;
end
c1 = c_grd_tm(1,end);
while c1 < nc_img
    c1 = c1+c_spc;
end

[c_ext, r_ext] = meshgrid(c0:c_spc:c1, r0:r_spc:r1);
[nr_ext, nc_ext] = size(c_ext);

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
u_ext_tx = spline2d(c_ext(:), r_ext(:), c_pts_tx(roi), r_pts_tx(roi), u_grd_tm(roi), tension);
v_ext_tx = spline2d(c_ext(:), r_ext(:), c_pts_tx(roi), r_pts_tx(roi), v_grd_tm(roi), tension);
u_ext_tx = reshape(u_ext_tx, nr_ext, nc_ext);
v_ext_tx = reshape(v_ext_tx, nr_ext, nc_ext);

% interpolate on low-res to full-res grid using cheap linear interpolation
interp = griddedInterpolant(r_ext, c_ext, u_ext_tx, 'linear');
u_img_tx = interp(r_img, c_img);
interp.Values = v_ext_tx;
v_img_tx = interp(r_img, c_img);

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