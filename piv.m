function [x_grd_tm, y_grd_tm, u_grd_tm, v_grd_tm, roi_grd_tm] = piv(...
    ini_ti, fin_tf, ini_roi_ti, fin_roi_tf, xw, yw, samplen, sampspc, ...
    intrlen, npass, valid_max, valid_eps, spline_tension, min_frac_data, ...
    min_frac_overlap, verbose)                 
% PIV analysis for Yalebox image data
%
% Arguments, input:
%
%   ini_ti, fin_tf = 2D matrix, double, range 0 to 1, normalized grayscale image from
%       the start and end of the step to be analyzed.
%
%   ini_roi_ti, fin_roi_tf = 2D matrix, logical, mask indicating pixels where there is
%       sand (1) and where there is only background (0) that should be ignored.
%
%   xw, yw = Vector, double, increasing, x- and y-direction coordinate vectors,
%       length must match the cols and rows, respectively, in both ini and fin.
%
%   samplen = Vector, length == number of grid resolutions, integer, side
%       length of the square sample window, [pixels]
%
%   sampspc = Scalar, integer, spacing between adjacent sample points in the
%       (square) sample grid, [pixels].
%
%   intrlen = Vector, length == number of grid resolutions, integer, side
%       length of the square interrogation window
%
%   npass = Vector, length == number of grid resolutions, integer, number of
%       image deformation passes
%
%   valid_max = Scalar, double, maximum value for the normalized residual
%       in the vector validation function, above which a vector is flagged
%       as invalid. Ref [3] reccomends a value of 2.
%
%   epsilon = Scalar, double, minumum value of the normalization factor in
%       the vector validation function. Ref [3] reccomends a value of 0.1.
%
%   spline_tension = Scalar, tension parameter for the spline interpolation
%       routine in ref [4]
%
%   min_frac_data = Scalar, minimum fraction of the sample window that must
%       contain data (e.g. sand) for the point to be included in the ROI for PIV
%       analysis
%
%   min_frac_overlap = Scalar, minimum fraction of the sample window data that
%       must overlap the interrogation window data for a point in the
%       cross-correlation to be valid
%
%   verbose = Scalar, integer, flag to enable (1) or diasable (0) verbose text
%       output messages
%
% Arguments, output:
%
%   xx, yy = Vector, double, coordinate vectors for the final output sample
%       grid, in world coordinate units
%
%   uu, vv = 2D matrix, double, computed displacement in the x- and y-directions
%       in world coordinate units
%
%   roi = 2D matrix, logical, flag indicating whether the data point lies within
%       PIV analysis (1) or not (0). Points inside the ROI may be measured or
%       interpolated, points outside the ROI are all interpolated/extrapolated.
%
% Notes:
%   + Sample grid spacing is held constant since our experiments showed
%     upsampling introduced ugly interpolation artifacts. We found that holding
%     grid spacing constant and reducing sample window size is an effective way
%     to increase spatial resolution while avoiding this issue.
%
% References:
%
% [1] Raffel, M., Willert, C. E., Wereley, S. T., & Kompenhans, J. (2007).
%   Particle Image Velocimetry: A Practical Guide. BOOK. 
%
% [2] Wereley, S. T. (2001). Adaptive Second-Order Accurate Particle Image
%   Velocimetry, 31
%
% [3] Westerweel, J., & Scarano, F. (2005). Universal outlier detection for PIV
%   data. Experiments in Fluids, 39(6), 1096???1100. doi:10.1007/s00348-005-0016-6
%
% [4] Wessel, P., & Bercovici, D. (1998). Interpolation with splines in tension:
%   A Green’s function approach. Mathematical Geology, 30(1), 77–93. Retrieved
%   from http://link.springer.com/article/10.1023/A:1021713421882

% EXPERIMENT toggle
% exp = true;
% exp = false;

% Note: variable suffixes are used to describe the time and space grids that
% each variable represents. These are:
%   ti -> time of the initial image
%   tm -> midpoint between initial and final images
%   tf -> time of the final image
%   grd -> regular sample grid
%   pts -> irregularly spaced points
%   img -> regular grid at image resolution
    
% check for sane inputs
[nr, nc] = size(ini_ti); % image size
ng = numel(samplen); % number of grid refinement steps
validateattributes( ini_ti,           {'double'},  {'2d', 'real', 'nonnan', '>=', 0, '<=' 1});
validateattributes( fin_tf,           {'double'},  {'2d', 'real', 'nonnan', '>=', 0, '<=' 1, 'size', [nr, nc]});
validateattributes( ini_roi_ti,       {'logical'}, {'2d', 'size', [nr, nc]});
validateattributes( fin_roi_tf,       {'logical'}, {'2d', 'size', [nr, nc]});
validateattributes( xw,               {'double'},  {'vector', 'real', 'nonnan', 'numel', nc});
validateattributes( yw,               {'double'},  {'vector', 'real', 'nonnan', 'numel', nr});
validateattributes( samplen,          {'numeric'}, {'vector', 'integer', 'positive', 'nonnan'});
validateattributes( sampspc,          {'numeric'}, {'scalar', 'integer', 'positive', 'nonnan'});
validateattributes( intrlen,          {'numeric'}, {'vector', 'numel', ng, 'integer', 'positive', 'nonnan'});
validateattributes( npass,            {'numeric'}, {'vector', 'numel', ng, 'integer', 'positive'});
validateattributes( valid_max,        {'double'},  {'scalar', 'positive'});
validateattributes( valid_eps,        {'double'},  {'scalar', 'positive'});
validateattributes( spline_tension,   {'numeric'}, {'scalar', '>=', 0, '<', 1});
validateattributes( min_frac_data,    {'numeric'}, {'scalar', '>=', 0, '<=', 1});
validateattributes( min_frac_overlap, {'numeric'}, {'scalar', '>=', 0, '<=', 1});
validateattributes( verbose,          {'numeric', 'logical'}, {'scalar', 'binary'});

% verbose output
if verbose
    fprintf('%s: ini: size = [%d, %d], roi frac = %.2f\n', mfilename, nr, nc, sum(ini_roi_ti(:))/numel(ini_roi_ti));
    fprintf('%s: fin: size = [%d, %d], roi frac = %.2f\n', mfilename, nr, nc, sum(fin_roi_tf(:))/numel(fin_roi_tf));
    fprintf('%s: xw: min = %.3f, max = %.3f\n', mfilename, min(xw), max(xw));
    fprintf('%s: yw: min = %.3f, max = %.3f\n', mfilename, min(yw), max(yw));
    fprintf('%s: samplen = ', mfilename); fprintf('%d ', samplen); fprintf('\n');
    fprintf('%s: sampspc = ', mfilename); fprintf('%d ', sampspc); fprintf('\n');
    fprintf('%s: intrlen = ', mfilename); fprintf('%d ', intrlen); fprintf('\n');
    fprintf('%s: npass = ', mfilename); fprintf('%d ', npass); fprintf('\n');
    fprintf('%s: valid_max = %.2f\n', mfilename, valid_max);
    fprintf('%s: valid_eps = %.2e\n', mfilename, valid_eps);
    fprintf('%s: spline_tension = %.3f\n', mfilename, spline_tension);
    fprintf('%s: min_frac_data = %.3f\n', mfilename, min_frac_data);
    fprintf('%s: min_frac_overlap = %.3f\n', mfilename, min_frac_overlap);
end

% TODO: it would make a LOT of sense to trim then restore the domain to its
% minimum size, this would reduce the cost of the full-image interpolations
% significantly

% expand grid definition vectors to reflect the number of passes
[samplen, intrlen] = expand_grid_def(samplen, intrlen, npass);

% init sample grid and solution
% NOTE: expect I will interpolate to this grid after the last step
[r_grd_tm, c_grd_tm, x_grd_tm, y_grd_tm] = piv_sample_grid(sampspc, xw, yw);
u_grd_tm = zeros(size(r_grd_tm));
v_grd_tm = zeros(size(c_grd_tm));

% multipass loop
np = length(samplen);
for pp = 1:np
    
    if verbose
        fprintf('%s: pass %d of %d: samplen = %d, intrlen = %d\n', ...
            mfilename, pp, np, samplen(pp), intrlen(pp));
    end
    
    % % jiggle sample points by random 1/4 sample spacing
    % r_pts = r_grd_tm + 0.50*sampspc*(rand(size(r_grd_tm)) - 0.5);
    % c_pts = c_grd_tm + 0.50*sampspc*(rand(size(c_grd_tm)) - 0.5);
    
    % get displacement update using normalized cross correlation
    if pp == np
        quality = true;
    else
        quality = false;
    end
    [r_pts_tm, c_pts_tm, u_pts_tm, v_pts_tm, roi_pts_tm] = piv_displacement_shift(...
        ini_ti, fin_tf, r_grd_tm, c_grd_tm, u_grd_tm, v_grd_tm, ...
        samplen(pp), intrlen(pp), min_frac_data, min_frac_overlap, quality, ...
        verbose);
    
    % validate displacement update
    % NOTE: neighborhood is hard-coded here
    % TODO: Neighborhood cutoff should use distance, some multiple of grid spacing
    % TODO: Needs a general review
    [u_pts_tm, v_pts_tm] = piv_validate_pts_nmed(...
        c_pts_tm, r_pts_tm, u_pts_tm, v_pts_tm, 8, valid_max, valid_eps, verbose);
     
%     % get roi for sample grid
%     % TODO: recompute ROI as alpha shape?
%     % TODO: see end, and delete if not used
%     roi_grd_tm = roi_pts_tm; % the lazy way

    % TODO: Cheap interpolation should be fine until the last pass - since the
    % solution is not accumulated

    % interpolate valid vectors to full sample grid (expensive)
    [u_grd_tm, v_grd_tm] = piv_interp_spline(...
        c_pts_tm, r_pts_tm, u_pts_tm, v_pts_tm, roi_pts_tm, ...
        c_grd_tm, r_grd_tm, true, spline_tension, verbose);
    
%     % interpolate valid vectors to the last ROI
%     % NOTE: outside ROI is set to 0, this avoids crazy guesses for later passes
%     % TODO: recompute ROI as alpha shape?
%     [u_grd_tm, v_grd_tm] = piv_interp_spline(...
%         c_pts_tm, r_pts_tm, u_pts_tm, v_pts_tm, roi_pts_tm, ...
%         c_grd_tm, r_grd_tm, roi_pts_tm, spline_tension, verbose);
    
end
% end multipass loop

% convert displacements to world coordinates (assumes equal grid spacing) 
u_grd_tm = u_grd_tm.*(xw(2) - xw(1));
v_grd_tm = v_grd_tm.*(yw(2) - yw(1));
roi_grd_tm = roi_pts_tm; % use the final ROI

% re-apply last ROI
% TODO: does this make sense? I think the strain routine just interpolates again.
u_grd_tm(~roi_grd_tm) = NaN;
v_grd_tm(~roi_grd_tm) = NaN;

end

%% subroutines

function [slen_ex, ilen_ex] = expand_grid_def(slen, ilen, np)
%
% Expand the grid definition vectors to include the correct number of passes for
% each grid. Input arguments are defined above, but use shortened names here:
% samplen -> slen, intrlen -> ilen, npass -> np. 
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