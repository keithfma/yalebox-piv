function [uu, vv, qual] = piv_validate_pts_nmed(...
    cc, rr, uu, vv, qual, radius, max_norm_res, epsilon)
%
% Validate the displacement vector field on a scattered grid using a
% normalized median test including all neighbors within a fixed radius.
% Returns a reduced set of only the valid points. See reference [1] for details
% on the original method for a regular grid.
%
% Arguments:
%   cc, rr = 2D matrix, double, location of displacement observations
%   uu, vv = 2D matrix, double, displacement vector components
%   qual: Matrix, Quality flags indicating how each observation was
%       computed (or not computed), possible values are enumerated in the
%       Quality class
%   radius = Scalar, distance in pixel units around each point to search
%       for neighbors to include in normalized median test
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
validateattributes(qual, {'Quality'}, {'size', size(uu)});
validateattributes(cc, {'numeric'}, {'size', size(uu)});
validateattributes(rr, {'numeric'}, {'size', size(uu)});

% initialize
num_nbr = nan(size(uu));
num_valid_init = sum(qual(:) == Quality.Valid);
kdtree = KDTreeSearcher([cc(:), rr(:)]);

% check all points
for kk = 1:numel(uu)
    
    % skip if there is no observation at this point
    if qual(kk) ~= Quality.Valid
        continue;
    end
    
    % get (valid) neighbors
    nbr = rangesearch(kdtree, [cc(kk), rr(kk)], radius);
    ind_nbr = cell2mat(nbr);
    ind_nbr(qual(ind_nbr) ~= Quality.Valid | ind_nbr == kk) = []; % drop self and invalid neighbors
    num_nbr(kk) = length(ind_nbr);
    unbr = uu(ind_nbr);
    vnbr = vv(ind_nbr);
    
    % compute neighbor median, residual, and median residual
    med_unbr = median(unbr);
    res_unbr = abs(unbr-med_unbr);
    med_res_unbr = median(res_unbr);
    
    med_vnbr = median(vnbr);
    res_vnbr = abs(vnbr-med_vnbr);
    med_res_vnbr = median(res_vnbr);
    
    % compute center normalized residual
    norm_res_u = abs(uu(kk)-med_unbr)/(med_res_unbr+epsilon);
    norm_res_v = abs(vv(kk)-med_vnbr)/(med_res_vnbr+epsilon);
    
    % combine vector components (max or sum)
    norm_res = max(norm_res_u, norm_res_v);
    
    % set quality flag if point is invalidated by this test
    if norm_res > max_norm_res
        qual(kk) = Quality.ValidationFailed;
    end
    
end

num_valid_fin = sum(qual(:) == Quality.Valid);
num_invalidated = num_valid_init - num_valid_fin;
pct_invalidated = 100*num_invalidated/num_valid_init;
fprintf('%s: invalidated %d (%.1f%%), mean num neighbors: %.2f\n', ...
    mfilename, num_invalidated, pct_invalidated, nanmean(num_nbr(:)));
