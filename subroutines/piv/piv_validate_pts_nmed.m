function [uu, vv] = piv_validate_pts_nmed(...
    cc, rr, uu, vv, radius, max_norm_res, epsilon)
%
% Validate the displacement vector field on a scattered grid using a
% normalized median test including all neighbors within a fixed radius.
% Returns a reduced set of only the valid points. See reference [1] for details
% on the original method for a regular grid.
%
% Arguments:
%   cc, rr = 2D matrix, double, location of scattered points for displacement
%       vectors.
%   uu, vv = 2D matrix, double, displacement vector components. NaNs are
%       treated as missing values to allow for roi masking.
%   radius = Scalar, distance around each point to search for neighbors to
%       include in normalized median test
%   max_norm_res = Scalar, double, maximum value for the normalized residual,
%       above which a vector is flagged as invalid. Reference [1] reccomends a
%       value of 2.
%   epsilon = Scalar, double, minumum value of the normalization factor.
%       Reference [1] recommends a value of 0.1.
%
% References:
%   [1] Westerweel, J., & Scarano, F. (2005). Universal outlier detection for
%   PIV data. Experiments in Fluids, 39(6), 1096???1100.
%   doi:10.1007/s00348-005-0016-6
% % 

% check inputs
validateattributes(uu, {'numeric'}, {});
validateattributes(vv, {'numeric'}, {'size', size(uu)});
validateattributes(cc, {'numeric'}, {'size', size(uu)});
validateattributes(rr, {'numeric'}, {'size', size(uu)});

% initialize
valid = true(size(uu));
num_nbr = nan(size(uu));
num_valid_init = sum(~isnan(uu(:)));
kdtree = KDTreeSearcher([cc(:), rr(:)]);

% check all points
% TODO: parfor yielded garbage results, not sure why
for kk = 1:numel(uu)
    
    % skip if ther e is no observation at this point
    if isnan(uu(kk))
        valid(kk) = false;
        continue;
    end
    
    % get neighbors
    nbr = rangesearch(kdtree, [cc(kk), rr(kk)], radius);
    ind_nbr = cell2mat(nbr);
    ind_nbr(ind_nbr == kk) = [];
    num_nbr(kk) = length(ind_nbr);
    unbr = uu(ind_nbr);
    vnbr = vv(ind_nbr);
    
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
    
    % mark invalid points, not set to NaN yet
    if norm_res > max_norm_res
        valid(kk) = false;
    end
    
end

% set invalid observations to NaN
uu(~valid) = NaN;
vv(~valid) = NaN;

num_valid_fin = sum(~isnan(uu(:)));
num_invalidated = num_valid_init - num_valid_fin;
pct_invalidated = 100*num_invalidated/num_valid_init;
fprintf('%s: invalidated %d (%.1f%%), mean num neighbors: %.2f\n', ...
    mfilename, num_invalidated, pct_invalidated, nanmean(num_nbr(:)));
