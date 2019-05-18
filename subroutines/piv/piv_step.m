function result = piv_step(...
    img_ti, img_tf, roi_img_ti, roi_img_tf, x_img, y_img, ...
    samp_len, samp_spc, intr_len, num_pass, valid_radius, valid_max, ...
    valid_eps, min_frac_data, min_frac_overlap)
% function result = piv_step(...
%     img_ti, img_tf, roi_img_ti, roi_img_tf, x_img, y_img, ...
%     samp_len, samp_spc, intr_len, num_pass, valid_radius, valid_max, ...
%     valid_eps, min_frac_data, min_frac_overlap)
%
% PIV analysis for Yalebox image data. Returns results at scattered (raw) and
% gridded (interpolated) points, evaluated at the midpoint time between the two
% input images. See subroutines docs for additional details.
%
% Arguments:
%   img_ti, img_tf = 2D matrix, double, range 0 to 1, normalized grayscale image
%       from the start and end of the step to be analyzed.
% 
%   roi_img_ti, roi_img_tf = 2D matrix, logical, mask indicating pixels where
%       there is sand (1) and background (0) that should be ignored.
%  
%   x_img, y_img = Vector, double, increasing, x- and y-direction
%       coordinate vectors, length must match the cols and rows, respectively,
%       in both image grids.
% 
%   samp_len = Vector, length == number of grid resolutions, integer, side
%       length of the square sample window, [pixels]
% 
%   samp_spc = Scalar, integer, spacing between adjacent sample points in the
%       (square) sample grid, [pixels].
% 
%   intr_len = Vector, length == number of grid resolutions, integer, side
%       length of the square interrogation window
% 
%   num_pass = Vector, length == number of grid resolutions, integer, number of
%       image deformation passes
% 
%   valid_radius: Scalar, radius around each sample point to include in vector
%       validation, recommended to use ~4*samp_spc, [pixel] units
% 
%   valid_max = Scalar, double, maximum value for the normalized residual
%       in the vector validation function, above which a vector is flagged
%       as invalid. Ref [3] reccomends a value of 2.
% 
%   valid_eps = Scalar, double, minimum value of the normalization factor in
%       the vector validation function. Ref [3] reccomends a value of 0.1.
% 
%   min_frac_data = Scalar, minimum fraction of the sample window that must
%       contain data (e.g. sand) for the point to be included in the ROI for PIV
%       analysis
% 
%   min_frac_overlap = Scalar, minimum fraction of the sample window data that
%       must overlap the interrogation window data for a point in the
%       cross-correlation to be valid
%
%  Returns:
% 
%   result.x, result.y: 2D matrix, x- and y-coordinates for gridded
%       sample points
% 
%   result.u, result.v: 2D matrix, x- and y-direction displacements
%       for gridded sample points, at initial time
% 
%   result.valid: 2D matrix, logical grid with set to "true" for
%       valid observations, and "false" for invalid
%
% Notes:
%   + Sample grid spacing is held constant since our experiments showed
%     upsampling introduced ugly interpolation artifacts. We found that
%     holding grid spacing constant and reducing sample window size is an
%     effective way to increase spatial resolution while avoiding this
%     issue.
%
% References:
%
% [1] Raffel, M., Willert, C. E., Wereley, S. T., & Kompenhans, J. (2007).
%   Particle Image Velocimetry: A Practical Guide. BOOK. 
% %

% NOTE: special suffixes describe the time and space grids, these are:
%   ti -> time of the initial image
%   tf -> time of the final image
%   grd -> regular sample grid
%   img -> regular grid at image resolution
% %

update_path('normxcorr2_masked', 'inpaint_nans');

% check for sane inputs
[nr, nc, nb] = size(img_ti); % image size
ng = numel(samp_len); % number of grid refinement steps
validateattributes(img_ti, {'uint8'},  {'3d'});
validateattributes(img_tf, {'uint8'},  {'3d', 'size', [nr, nc, nb]});
validateattributes(roi_img_ti, {'logical'}, {'2d', 'size', [nr, nc]});
validateattributes(roi_img_tf, {'logical'}, {'2d', 'size', [nr, nc]});
validateattributes(x_img, {'double'},  {'vector', 'real', 'nonnan', 'numel', nc});
validateattributes(y_img, {'double'},  {'vector', 'real', 'nonnan', 'numel', nr});
validateattributes(samp_len, {'numeric'}, {'vector', 'integer', 'positive', 'nonnan'});
validateattributes(samp_spc, {'numeric'}, {'scalar', 'integer', 'positive', 'nonnan'});
validateattributes(intr_len, {'numeric'}, {'vector', 'numel', ng, 'integer', 'positive', 'nonnan'});
validateattributes(num_pass, {'numeric'}, {'vector', 'numel', ng, 'integer', 'positive'});
validateattributes(valid_radius, {'numeric'}, {'scalar', 'integer', 'positive', 'nonnan'});
validateattributes(valid_max, {'double'}, {'scalar', 'positive'});
validateattributes(valid_eps, {'double'}, {'scalar', 'positive'});
validateattributes(min_frac_data, {'numeric'}, {'scalar', '>=', 0, '<=', 1});
validateattributes(min_frac_overlap, {'numeric'}, {'scalar', '>=', 0, '<=', 1});

fprintf('%s: ini: size = [%d, %d], roi frac = %.2f\n', mfilename, nr, nc, sum(roi_img_ti(:))/numel(roi_img_ti));
fprintf('%s: fin: size = [%d, %d], roi frac = %.2f\n', mfilename, nr, nc, sum(roi_img_tf(:))/numel(roi_img_tf));
fprintf('%s: x_img: min = %.3f, max = %.3f\n', mfilename, min(x_img), max(x_img));
fprintf('%s: y_img: min = %.3f, max = %.3f\n', mfilename, min(y_img), max(y_img));
fprintf('%s: samp_len = ', mfilename); fprintf('%d ', samp_len); fprintf('\n');
fprintf('%s: samp_spc = ', mfilename); fprintf('%d ', samp_spc); fprintf('\n');
fprintf('%s: intr_len = ', mfilename); fprintf('%d ', intr_len); fprintf('\n');
fprintf('%s: num_pass = ', mfilename); fprintf('%d ', num_pass); fprintf('\n');
fprintf('%s: valid_radius = %.2e\n', mfilename, valid_radius);
fprintf('%s: valid_max = %.2f\n', mfilename, valid_max);
fprintf('%s: valid_eps = %.2e\n', mfilename, valid_eps);
fprintf('%s: min_frac_data = %.3f\n', mfilename, min_frac_data);
fprintf('%s: min_frac_overlap = %.3f\n', mfilename, min_frac_overlap);

% treat image boundaries
% TODO: choose one, and drop other options
boundary_treatment = 'extend';
fprintf('%s: apply boundary treatment "%s" to input images\n', ...
    mfilename, boundary_treatment); 
[~, ~, img_ti] = piv_fill_nnbr(...
    x_img, y_img, img_ti, roi_img_ti, samp_len(1), boundary_treatment);
[x_img, y_img, img_tf] = piv_fill_nnbr(...
    x_img, y_img, img_tf, roi_img_ti, samp_len(1), boundary_treatment);

% convert to images to grayscale and apply mask
img_ti = rgb2hsv(img_ti);
img_ti = img_ti(:,:,3);

img_tf = rgb2hsv(img_tf);
img_tf = img_tf(:,:,3);

% DEBUG: save a gif to the current working directory
gif_file = sprintf('test_%s.gif', boundary_treatment);
imwrite(uint8(img_ti*255), gif_file, 'gif', 'Loopcount', inf, 'DelayTime', 1);
imwrite(uint8(img_tf*255), gif_file, 'gif','WriteMode','append', 'DelayTime', 1);

% create sample grid
[rr, cc] = piv_sample_grid(samp_len, samp_spc, size(img_ti, 1), size(img_ti, 2));

% expand grid definition vectors to reflect the number of passes
[samp_len, intr_len] = expand_grid_def(samp_len, intr_len, num_pass);

% initial guess for displacements
uu = zeros(size(rr));
vv = zeros(size(cc));
    
% multipass loop
np = length(samp_len);
for pp = 1:np
    
    fprintf('%s: pass %d of %d: samp_len = %d, intr_len = %d\n', ...
        mfilename, pp, np, samp_len(pp), intr_len(pp));
    
    % get displacement update using normalized cross correlation
    [uu, vv] = piv_displacement(...
        img_ti, img_tf, rr, cc, uu, vv, samp_len(pp), ...
        intr_len(pp), min_frac_data, min_frac_overlap);
    
    % validate displacement vectors
    [uu, vv] = piv_validate_pts_nmed(...
        cc, rr, uu, vv, valid_radius, valid_max, valid_eps);
     
    % interpolate valid vectors to full sample grid 
    % note: leave as-is on last pass, want NaN where there is no valid
    %   observation, these can be filled later if desired, but should not
    %   pretend to be good data
    if pp < np
        uu = inpaint_nans(uu, 2);
        vv = inpaint_nans(vv, 2);
    end
    
end
% end multipass loop

% convert all output variables to world coordinate system
[xx, yy] = coord_intrinsic_to_world(rr, cc, x_img, y_img);
[uu, vv] = displ_intrinsic_to_world(uu, vv, x_img, y_img);

% package results as structure
result = struct(...
    'x', xx, 'y', yy, 'u', uu, 'v', vv);

end

%% subroutines

function [xx, yy] = coord_intrinsic_to_world(rr, cc, xx_ref, yy_ref)
%
% Convert intrinsic (pixel) coordinates to world coordinates
%
% Arguments:
%   rr, cc: Row and column intrinsic coordinates to be converted
%   xx_ref, yy_ref: Vectors, reference coordinate vectors
%   xx, yy: World coordinates correspondint to rr, cc
%
% Note: Assumes that xx_ref, yy_ref define a regular coordinate grid with no
%   rotation (i.e. x = f(r, c) = f(c)), and that they are in the same intrinsic
%   coordinate system as rr, cc (i.e. the world coordinate position of the point
%   rr = 1, cc = 1 is exactly xx = xx_ref(1), yy = yy_ref(1))
% % 

xx = interp1(1:length(xx_ref), xx_ref, cc, 'linear', 'extrap');
yy = interp1(1:length(yy_ref), yy_ref, rr, 'linear', 'extrap');

end

function [uw, vw] = displ_intrinsic_to_world(ui, vi, xx_ref, yy_ref)
%
% Convert displacements from intrisic (pixels) to world units
%
% Arguments:
%   ui, vi: Row and column displacements in intrinsic coordinates to convert
%   xx_ref, yy_ref: Vectors, reference coordinate vectors
%   uf, vf: x- and y- direction displacements in world coordinates
%
% Note: Assumes that xx_ref, yy_ref define a regular coordinate grid with no
%   rotation.
% % 

dx = (xx_ref(2) - xx_ref(1));
uw = ui*dx;

dy = (yy_ref(2) - yy_ref(1));
vw = vi*dy;

end

function [slen_ex, ilen_ex] = expand_grid_def(slen, ilen, np)
%
% Expand the grid definition vectors to include the correct number of passes for
% each grid. Input arguments are defined above, but use shortened names here:
% samp_len -> slen, intr_len -> ilen, num_pass -> np. 
%
% Note: outputs are intentionally not preallocated - these vectors are small and
% the performace cost is negligible. 
% %

slen_ex = [];
ilen_ex = [];

for ii = 1:length(np)
   slen_ex = [slen_ex, repmat(slen(ii), 1, np(ii))]; %#ok!
   ilen_ex = [ilen_ex, repmat(ilen(ii), 1, np(ii))]; %#ok!
end

end
