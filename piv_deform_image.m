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

% NOTE: roi is currently unused

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

% interpolate vectors to full image resolution as complex field (cheaply)
[u_img_tx, v_img_tx] = piv_interp_linear(c_pts_tx, r_pts_tx, ...
    u_grd_tm, v_grd_tm, c_img, r_img, true, verbose);

% prep displacement matrix for image deformation 
if is_fwd
    displacement = -0.5*cat(3, u_img_tx, v_img_tx);
else
    displacement = 0.5*cat(3, u_img_tx, v_img_tx);
end

% deform image to midpoint time
img_tm = imwarp(img_tx, displacement, 'cubic', 'FillValues', 0);

% deform roi to midpoint time
% TODO: add ROI as another layer, and warp once.
tmp = imwarp(double(img_roi_tx), displacement, 'cubic', 'FillValues', 0);
img_roi_tm = abs(tmp-1) < roi_epsilon;

% re-apply roi
img_tm(~img_roi_tm) = 0;