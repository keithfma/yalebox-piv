function result = piv_step(...
    img_ti, img_tf, mask_ti, mask_tf, x_img, y_img, ...
    samp_len, samp_spc, intr_len, num_pass, valid_radius, valid_max, ...
    valid_eps, min_frac_data, min_frac_overlap)
% function result = piv_step(...
%     img_ti, img_tf, mask_ti, mask_tf, x_img, y_img, ...
%     samp_len, samp_spc, intr_len, num_pass, valid_radius, valid_max, ...
%     valid_eps, min_frac_data, min_frac_overlap)
%
% PIV analysis for Yalebox image data. Returns results evaluated at the
% time of the initial image. See subroutines docs for additional details.
%
% Arguments:
%   img_ti, img_tf = 2D matrix, double, range 0 to 1, normalized grayscale
%       images from the start and end of the step to be analyzed.
%   mask_ti, mask_tf = 2D matrix, logical, mask indicating pixels where
%       there is sand (1) and background (0) that should be ignored.
%   x_img, y_img = Vector, double, increasing, x- and y-direction
%       coordinate vectors, length must match the cols and rows, respectively,
%       in both image grids.
%   samp_len = Vector, length == number of grid resolutions, integer, side
%       length of the square sample window, [pixels]
%   samp_spc = Scalar, integer, spacing between adjacent sample points in the
%       (square) sample grid, [pixels].
%   intr_len = Vector, length == number of grid resolutions, integer, side
%       length of the square interrogation window
%   num_pass = Vector, length == number of grid resolutions, integer, number of
%       image deformation passes
%   valid_radius: Scalar, radius around each sample point to include in vector
%       validation, recommended to use ~4*samp_spc, [pixel] units
%   valid_max = Scalar, double, maximum value for the normalized residual
%       in the vector validation function, above which a vector is flagged
%       as invalid. Ref [3] reccomends a value of 2.
%   valid_eps = Scalar, double, minimum value of the normalization factor in
%       the vector validation function. Ref [3] reccomends a value of 0.1.
%   min_frac_data = Scalar, minimum fraction of the sample window that must
%       contain data (e.g. sand) for the point to be included in the ROI for PIV
%       analysis 
%   min_frac_overlap = Scalar, minimum fraction of the sample window data that
%       must overlap the interrogation window data for a point in the
%       cross-correlation to be valid
%
%  Returns:
% 
%   result.x, result.y: 2D matrix, x- and y-coordinates for gridded
%       sample points
%   result.u, result.v: 2D matrix, x- and y-direction displacements
%       for gridded sample points, at initial time
%   result.quality: 2D matrix, integer labels indicating how this
%       observation was computed, possible values are enumerated in the
%       Quality class
%   result.mask: 2D matrix, logical labels indicating if an observation is
%       inside (true) or outside (false) the original sand region, a
%       given observation is considered "inside" if a majority of its
%       sample window pixels are within the image mask.
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
%   img -> regular grid at image resolution
% %

update_path('normxcorr2_masked', 'inpaint_nans');

% check for sane inputs
[nr, nc] = size(img_ti); % image size
ng = numel(samp_len); % number of grid refinement steps
validateattributes(img_ti, {'double', 'single'},  {'2d'});
validateattributes(img_tf, {'double', 'single'},  {'2d', 'size', [nr, nc]});
validateattributes(mask_ti, {'logical'}, {'2d', 'size', [nr, nc]});
validateattributes(mask_tf, {'logical'}, {'2d', 'size', [nr, nc]});
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

fprintf('%s: ini: size = [%d, %d], mask frac = %.2f\n', mfilename, nr, nc, sum(mask_ti(:))/numel(mask_ti));
fprintf('%s: fin: size = [%d, %d], mask frac = %.2f\n', mfilename, nr, nc, sum(mask_tf(:))/numel(mask_tf));
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

% create sample grid and pad image arrays
% FIXME: update piv_sample_grid such that the grid starts at the origin, this
%   change is needed to respect the fact that fill has been done elsewhere
[rr, cc] = piv_sample_grid(samp_len, samp_spc, x_img, y_img);

% expand grid definition vectors to reflect the number of passes
[samp_len, intr_len] = expand_grid_def(samp_len, intr_len, num_pass);

% create mask array by interpolating the image masks down to the sample grid,
%   initial points and final points must be within the image 
%   mask to be labled as true
[cc_img, rr_img] = meshgrid(1:size(img_ti, 2), 1:size(img_ti, 1));
samp_mask_ti = interp2(cc_img, rr_img, single(mask_ti), cc, rr, 'nearest'); 
samp_mask_tf = interp2(cc_img, rr_img, single(mask_tf), cc, rr, 'nearest'); 
mask = samp_mask_ti & samp_mask_tf;

% initial guess for displacements
uu = zeros(size(rr));
vv = zeros(size(cc));
    
% multipass loop
np = length(samp_len);
for pp = 1:np
    
    fprintf('%s: pass %d of %d: samp_len = %d, intr_len = %d\n', ...
        mfilename, pp, np, samp_len(pp), intr_len(pp));
    
    % get displacement update using normalized cross correlation
    [uu, vv, qual] = piv_displacement(...
        img_ti, img_tf, rr, cc, uu, vv, samp_len(pp), ...
        intr_len(pp), min_frac_data, min_frac_overlap);
    
    if pp == 1
        % drop fill region for initial pass
        % note: reason is to populate a decent initial guess by interpolation
        qual(~mask) = Quality.Skipped;
    end
    
    % validate displacement vectors
    [uu, vv, qual] = piv_validate_pts_nmed(...
        cc, rr, uu, vv, qual, valid_radius, valid_max, valid_eps);
     
    % interpolate valid vectors to full sample grid
    valid = qual == Quality.Valid;
    interpolant = scatteredInterpolant(cc(valid), rr(valid), uu(valid), 'natural', 'nearest');
    uu(~valid) = interpolant(cc(~valid), rr(~valid));
    interpolant.Values = vv(valid);
    vv(~valid) = interpolant(cc(~valid), rr(~valid));
    
end
% end multipass loop

% convert all output variables to world coordinate system
[xx, yy] = coord_intrinsic_to_world(rr, cc, x_img, y_img);
[uu, vv] = displ_intrinsic_to_world(uu, vv, x_img, y_img);

% package results as structure
x = xx(1, :); x = x(:);  % convert to vectors
y = yy(:, 1); y = y(:);
result = struct(...
    'x', x, 'y', y, 'u', uu, 'v', vv, 'quality', qual, 'mask', mask);

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
