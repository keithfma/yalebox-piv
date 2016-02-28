function [uu, vv] = piv_validate_pts_nmed(cc, rr, uu, vv, num_nbr, max_norm_res, epsilon, verbose)
%
% Validate the displacement vector field on a scattered grid using a normalized
% median test including the nearest num_nbr neighbors.  Invalidated points are
% set to NaN. See reference [1] for details on the original method for a regular
% grid.
%
% Arguments:
%
%   cc, rr = 2D matrix, double, location of scattered points for displacement
%       vectors.
%
%   uu, vv = 2D matrix, double, displacement vector components. NaNs are
%       treated as missing values to allow for roi masking.
%
%   num_nbr = Scalar integer, number of nearest neighbors to include in
%       normalized median test
%
%   max_norm_res = Scalar, double, maximum value for the normalized residual,
%       above which a vector is flagged as invalid. Reference [1] reccomends a
%       value of 2.
%
%   epsilon = Scalar, double, minumum value of the normalization factor.
%       Reference [1] recommends a value of 0.1.
%
%   verbose = Enable (1) or disable (2) verbose output
%
% References:
%
% [1] Westerweel, J., & Scarano, F. (2005). Universal outlier detection for PIV
% data. Experiments in Fluids, 39(6), 1096???1100. doi:10.1007/s00348-005-0016-6

% initialize
valid = ~isnan(vv) & ~isnan(uu);
num_valid_init = sum(valid(:));
kdtree = KDTreeSearcher([cc(valid), rr(valid)]);

ind = reshape(1:numel(uu), size(uu));
ind = ind(valid)'; % linear index for initially valid points, must be a row vector, so we can loop over it

% loop over initially valid points only
for kk = ind
    
    % get neighbors
    nbr = knnsearch(kdtree, [cc(kk), rr(kk)], 'K', num_nbr+1, 'IncludeTies', true);
    ind_nbr = ind(cell2mat(nbr));
    ind_nbr(ind_nbr == kk) = [];
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

% set invalid points to NaN
uu(~valid) = NaN;
vv(~valid) = NaN;

if verbose
    num_invalidated = num_valid_init-sum(valid(:));
    pct_invalidated = 100*num_invalidated/num_valid_init;
    fprintf('%s: %d invalidated (%.1f%%)\n', mfilename, num_invalidated, pct_invalidated);
end
   