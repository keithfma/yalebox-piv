function [xx, yy, uu, vv, roi] = ...
    piv(ini_ti, fin_tf, ini_roi_ti, fin_roi_tf, xw, yw, samplen, sampspc, ...
        intrlen, npass, valid_max, valid_eps, lowess_span_pts, verbose)                 
% New implementation PIV analysis for Yalebox image data
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
%   lowess_span_pts = Scalar, number of points to include in local fit for 
%       lowess smoothing/interpolation step.
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

% Note: variable suffixes are used to describe the time and space grids that
% each variable represents. These are:
%   ti -> time of the initial image
%   tm -> midpoint between initial and final images
%   tf -> time of the final image
%   grd -> regular sample grid
%   pts -> irregularly spaced points
%   img -> regular grid at image resolution

% local parameters
tension = 0.95;
        
% parse inputs
check_input(ini_ti, fin_tf, ini_roi_ti, fin_roi_tf, xw, yw, samplen, sampspc, intrlen, ...
    npass, valid_max, valid_eps, lowess_span_pts, verbose);

% expand grid definition vectors to reflect the number of passes
[samplen, intrlen, sampspc] = expand_grid_def(samplen, intrlen, sampspc, npass);

% init coordinate grids
[r_grd, c_grd, xx, yy] = piv_sample_grid(samplen(1), sampspc(1), xw, yw);

% preallocate as reccomended by Mlint
sz = size(c_grd);
u_grd_tm = zeros(sz); 
v_grd_tm = zeros(sz);  

% init deformed images
ini_tm = ini_ti;
fin_tm = fin_tf;

% multipass loop
np = length(samplen); 
for pp = 1:np
    
    % get displacements update using normalized cross correlation
    [r_pts, c_pts, du_pts_tm, dv_pts_tm, roi] = ...
        piv_displacement(ini_tm, fin_tm, r_grd, c_grd, samplen(pp), intrlen(pp));
    
    % validate displacement update 
    [du_pts_tm, dv_pts_tm] = ...
        piv_validate_pts_nmed(c_pts, r_pts, du_pts_tm, dv_pts_tm, 8, valid_max, ...
            valid_eps, true);    

    % interpolate/smooth valid vectors to sample grid, outside roi is NaN
    [du_grd_tm, dv_grd_tm] = ...
        piv_lowess_interp(c_pts, r_pts, du_pts_tm, dv_pts_tm, c_grd, r_grd, roi, lowess_span_pts);
    
    % update displacement, points outside roi become NaN
    u_grd_tm = u_grd_tm + du_grd_tm;
    v_grd_tm = v_grd_tm + dv_grd_tm;
    
    % deform images to midpoint time, if there is another pass
    if pp < np
        ini_tm = piv_deform_image(ini_ti, ini_roi_ti, r_grd, c_grd, u_grd_tm, ...
            v_grd_tm, roi, 1);
        fin_tm = piv_deform_image(fin_tf, fin_roi_tf, r_grd, c_grd, u_grd_tm, ...
            v_grd_tm, roi, 0);        
    end
    
    % interpolate to new sample grid, if grid is changed in the next pass
    % ...interpolation fills the whole sample grid, the new roi will be imposed 
    % ...by adding NaNs to the current estimate 
    if pp<np && (samplen(pp)~=samplen(pp+1) || sampspc(pp)~=sampspc(pp+1))
        
        [r_grd_next, c_grd_next, xx, yy] = piv_sample_grid(samplen(pp+1), sampspc(pp+1), xw, yw);        
        u_grd_tm = spline2d(c_grd_next(:), r_grd_next(:), c_grd(roi), r_grd(roi), u_grd_tm(roi), tension);        
        v_grd_tm = spline2d(c_grd_next(:), r_grd_next(:), c_grd(roi), r_grd(roi), v_grd_tm(roi), tension);
        r_grd = r_grd_next;
        c_grd = c_grd_next;        
        u_grd_tm = reshape(u_grd_tm, size(r_grd));
        v_grd_tm = reshape(v_grd_tm, size(r_grd));
        
    end
    
end   
% end multipass loop

% convert displacements to world coordinates (assumes equal grid spacing)
uu = u_grd_tm.*(xw(2)-xw(1)); 
vv = v_grd_tm.*(yw(2)-yw(1));

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

function [] = check_input(ini, fin, ini_roi, fin_roi, xx, yy, samplen, ...
    sampspc, intrlen, npass, valid_max, valid_eps, lowess_span_pts, verbose)
%
% Check for sane input argument properties, exit with error if they do not
% match expectations.
% %

[nr, nc] = size(ini); % image size
ng = numel(samplen); % number of grid refinement steps

validateattributes(ini, {'double'}, {'2d', 'real', 'nonnan', '>=', 0, '<=' 1});
validateattributes(fin, {'double'}, {'2d', 'real', 'nonnan', '>=', 0, '<=' 1, ...
    'size', [nr, nc]});
validateattributes(ini_roi, {'logical'}, {'2d', 'size', [nr, nc]});
validateattributes(fin_roi, {'logical'}, {'2d', 'size', [nr, nc]});
validateattributes(xx, {'double'}, {'vector', 'real', 'nonnan', 'numel', nc});
validateattributes(yy, {'double'}, {'vector', 'real', 'nonnan', 'numel', nr});
validateattributes(samplen, {'numeric'}, {'vector', 'integer', 'positive', ...
    'nonnan'});
validateattributes(sampspc, {'numeric'}, {'vector', 'numel', ng, 'integer', ...
    'positive', 'nonnan'});
validateattributes(intrlen, {'numeric'}, {'vector', 'numel', ng, 'integer', ...
    'positive', 'nonnan'});
validateattributes(npass, {'numeric'}, {'vector', 'numel', ng, 'integer', ...
    'positive'});
validateattributes(valid_max, {'double'}, {'scalar', 'positive'});
validateattributes(valid_eps, {'double'}, {'scalar', 'positive'});
validateattributes(lowess_span_pts, {'numeric'}, {'scalar', 'integer'});
validateattributes(verbose, {'numeric', 'logical'}, {'scalar', 'binary'});

end
