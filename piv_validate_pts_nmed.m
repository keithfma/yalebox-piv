function [uu, vv] = piv_validate_pts_nmed(cc, rr, uu, vv, radius, max_norm_res, epsilon)
%
% Validate the displacement vector field using a normalized median test
% including all neighbors within a radius nbr_radius.  Invalidated points are
% set to NaN. See reference [1] for details.
%
% Arguments:
%
%   cc = 
%
%   rr = 
% 
%   uu, vv = 2D matrix, double, displacement vector components. NaNs are
%       treated as missing values to allow for roi masking.
%
%   radius = Scalar double, radius of circular neighborhood about each point
%       to be included in the test, default value sqrt(2) 
%
%   max_norm_res = Scalar, double, maximum value for the normalized residual,
%       above which a vector is flagged as invalid. Reference [1] reccomends a
%       value of 2.
%
%   epsilon = Scalar, double, minumum value of the normalization factor.
%       Reference [1] recommends a value of 0.1.
%
% References: 
%
% [1] Westerweel, J., & Scarano, F. (2005). Universal outlier detection for PIV
% data. Experiments in Fluids, 39(6), 1096???1100. doi:10.1007/s00348-005-0016-6

% Debug parameters
min_nbrs = 8;

% % set defaults
% if nargin < 3
%     nbr_radius = sqrt(2);
% end
% if nargin < 4; 
%     max_norm_res = 2; 
% end
% if nargin < 5; 
%     epsilon = 0.1; 
% end 

% precompute constants
outer_box_dist = radius;
inner_box_dist = radius*sqrt(2)/2;
radius_squared = radius^2;

% % drop all NaNs
% keep = ~isnan(uu) & ~isnan(vv);
% uu = uu(keep);
% vv = vv(keep);
% cc = cc(keep);
% rr = rr(keep);

% loop over all points
for kk = 1:numel(uu)
    
    % skip if point is already invalidated
    if isnan(uu(kk)) || isnan(vv(kk))
        continue
    end
 
    % get neighbors
    dc = abs(cc-cc(kk));
    dr = abs(rr-rr(kk));
    outer = (abs(dc) <= outer_box_dist) & (abs(dr) <= outer_box_dist);
    inner = (abs(dc) <= inner_box_dist) & (abs(dr) <= inner_box_dist);
    between = outer & ~inner;
    between( (dr(between).^2+dc(between).^2) > radius_squared ) = false;
    is_nbr = inner | between;    
    is_nbr(kk) = false;
    unbr = uu(is_nbr);
    vnbr = vv(is_nbr);
 
    % invalidate if point has too few neighbors (should not happen for smooth fields)
    fprintf('num_nbrs = %i\n', numel(unbr));
    
%     figure(1); 
%     hold off; 
%     plot(cc, rr, '.k');  
%     hold on; 
%     plot(cc(is_nbr), rr(is_nbr), '.b');    
%     plot(cc(kk), rr(kk), 'xr');    
    
    % NOTE: the count is incorrect since it includes NaNs
    
    if numel(unbr) < min_nbrs
        uu(kk) = NaN;
        vv(kk) = NaN;
        continue
    end

    % compute neighbor median, residual, and median residual
    med_unbr = nanmedian(unbr);
    res_unbr = abs(unbr-med_unbr);
    med_res_unbr = nanmedian(res_unbr);
    
    med_vnbr = nanmedian(vnbr);
    res_vnbr = abs(vnbr-med_vnbr);
    med_res_vnbr = nanmedian(res_vnbr);
    
    % compute center normalized residual
    norm_res_u = abs(uu(kk)-med_unbr)/(med_res_unbr+epsilon);
    norm_res_v = abs(vv(kk)-med_vnbr)/(med_res_vnbr+epsilon);
    
    % combine vector components (max or sum)
    norm_res = max(norm_res_u, norm_res_v);
    
    % set invalid points to NaN
    fprintf('norm_res = %f\n', norm_res);
    if norm_res > max_norm_res
        uu(kk) = NaN;
        vv(kk) = NaN;
    end
    
    % NOTE: the available neighbors are changing as they are set to NAN
    % NOTE: a multipass approach might be appropriate, iterate until no more are invalidated
    
end