function [x_soln, y_soln, u_soln_tm, v_soln_tm, roi_soln] = piv(...
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
%       length of the square sample window
%
%   sampspc = Vector, length == number of grid resolutions, integer,
%       spacing between adjacent sample points in the (square) sample grid
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
validateattributes( sampspc,          {'numeric'}, {'vector', 'numel', ng, 'integer', 'positive', 'nonnan'});
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

% NOTE: it would make a LOT of sense to trim then restore the domain to its
% minimum size, this would reduce the cost of the full-image interpolations
% significantly

% expand grid definition vectors to reflect the number of passes
[samplen, intrlen, sampspc] = expand_grid_def(samplen, intrlen, sampspc, npass);

% init coordinate grids for accumulated solution at final pass resolution
[r_soln, c_soln, x_soln, y_soln] = piv_sample_grid(...
    samplen(end), sampspc(end), xw, yw);
u_soln_tm = zeros(size(r_soln)); 
v_soln_tm = zeros(size(r_soln));  

% init deformed images
ini_tm = ini_ti;
fin_tm = fin_tf;

% multipass loop
np = length(samplen);
for pp = 1:np
    
    if verbose
        fprintf('%s: pass %d of %d: samplen = %d, sampspc = %d, intrlen = %d\n', ...
            mfilename, pp, np, samplen(pp), sampspc(pp), intrlen(pp));
    end

    % create sampling coordinate grids 
    [r_grd, c_grd, ~, ~] = piv_sample_grid(samplen(pp), sampspc(pp), xw, yw);
    
    % jiggle sample points by random 1/4 sample spacing
    r_smp = r_grd + 0.50*sampspc(pp)*(rand(size(r_grd)) - 0.5);
    c_smp = c_grd + 0.50*sampspc(pp)*(rand(size(c_grd)) - 0.5);
    
    % get displacement update using normalized cross correlation
    [r_pts, c_pts, du_pts_tm, dv_pts_tm, roi] = piv_displacement(...
        ini_tm, fin_tm, r_smp, c_smp, samplen(pp), intrlen(pp), ...
        min_frac_data, min_frac_overlap, verbose);
    
    % validate displacement update
    % NOTE: neighborhood is hard-coded here
    % NOTE: 0.1 seems too high for valid_eps
    % NOTE: 8 seems too low for # neighbors
    [du_pts_tm, dv_pts_tm] = piv_validate_pts_nmed(...
        c_pts, r_pts, du_pts_tm, dv_pts_tm, 16, valid_max, 0.05, verbose);
    
    % NOTE: can probably figure out a way to interpolate to only the ROI at
    %   higher resolution, or, could do a cheap within-alpha-shape interpolation
    
    % NOTE: this interpolation still introduces terrifying ringing
    
    % <EXPERIMENT>
    % NOTE: alpha is not set, probably should be some multiple of the grid spacing
    [du_soln_tm, dv_soln_tm] = piv_interp_linear_alpha(...
        c_pts, r_pts, du_pts_tm, dv_pts_tm, c_soln, r_soln, true, ...
        3*sampspc(pp), verbose);
    % </EXPERIMENT>
    
    % % <ORIGINAL>
    % % interpolate valid vectors to full accumulated solution grid (expensive)
    % [du_soln_tm, dv_soln_tm] = piv_interp_spline(...
    %     c_pts, r_pts, du_pts_tm, dv_pts_tm, c_soln, r_soln, true, ...
    %     spline_tension, verbose);
    % % <ORIGINAL>
    
    % update displacement solution
    u_soln_tm = u_soln_tm + du_soln_tm;
    v_soln_tm = v_soln_tm + dv_soln_tm;
    
    keyboard
    
    % prepare for next pass, if needed
    if pp < np
        
        % deform images to midpoint time
        % NOTE: roi variable is unused, set to []
        ini_tm = piv_deform_image(ini_ti, ini_roi_ti, r_soln, c_soln, ...
            u_soln_tm, v_soln_tm, [], spline_tension, 1, verbose);
        fin_tm = piv_deform_image(fin_tf, fin_roi_tf, r_soln, c_soln, ...
            u_soln_tm, v_soln_tm, [], spline_tension, 0, verbose); 
    end
    
end
% end multipass loop

% convert displacements to world coordinates (assumes equal grid spacing)
u_soln_tm = u_soln_tm.*(xw(2) - xw(1));
v_soln_tm = v_soln_tm.*(yw(2) - yw(1));
roi_soln = roi;

end

%% subroutines

function [slen_ex, ilen_ex, sspc_ex] = expand_grid_def(slen, ilen, sspc, np)
%
% Expand the grid definition vectors to include the correct number of passes for
% each grid. Input arguments are defined above, but use shortened names here:
% samplen -> slen, intrlen -> ilen, sampspc -> sspc, npass -> np. 
%
% Note: outputs are intentionally not preallocated - these vectors are small and
% the performace cost is negligible. 
% %

slen_ex = [];
ilen_ex = [];
sspc_ex = [];

for ii = 1:length(np)
   slen_ex = [slen_ex, repmat(slen(ii), 1, np(ii))]; %#ok!
   ilen_ex = [ilen_ex, repmat(ilen(ii), 1, np(ii))]; %#ok!
   sspc_ex = [sspc_ex, repmat(sspc(ii), 1, np(ii))]; %#ok!
end

end