function [uu, vv] = yalebox_piv_valid_nmed(uu, vv, num_nbr, max_norm_res, epsilon)
% function [uu, vv] = yalebox_piv_valid_nmed(uu, vv, max_norm_res, epsilon)
%
% Validate the displacement vector field using a normalized median test with a
% num_nbr-point neighborhood. Neighborhood includes nearest num_nbr non-NaN
% points to account for domain edges. Invalidated points are set to NaN. See
% reference [1] for details.
%
% Arguments:
%
%   uu, vv = 2D matrix, double, displacement vector components. NaNs are
%       treated as missing values to allow for roi masking.
%
%   num_nbr = Scalar integer, number of nearest non-NaN neighbors to include in
%       test, default value is 8
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

% define parameters
default_num_nbr = 8;
default_max_norm_res = 2;
default_epsilon = 0.1; 

% set defaults
if nargin < 3
    num_nbr = default_num_nbr;
end
if nargin < 4; 
    max_norm_res = default_max_norm_res; 
end
if nargin < 5; 
    epsilon = default_epsilon; 
end 

% init
[nr, nc] = size(uu);

% get offsets to neighbors sorted by distance to the central point
%... 7x7 neighborhood without central point
roffset = [-3; -2; -1; 0; 1; 2; 3]*[ 1,  1,  1, 1, 1, 1, 1];
coffset = [ 1;  1;  1; 1; 1; 1; 1]*[-3, -2, -1, 0, 1, 2, 3];
%... remove central point
keep_ind = [1:24, 26:49];
roffset = roffset(keep_ind);
coffset = coffset(keep_ind);
%... sort by distance to center (0,0)
doffset = sqrt(roffset.^2+coffset.^2);
[~, sort_ind] = sort(doffset);
roffset = roffset(sort_ind);
coffset = coffset(sort_ind);

% loop over all displacement vectors
for ii = 1:nr
    for jj = 1:nc
        
        % skip if center point is already invalidated
        if isnan(uu(ii,jj)) || isnan(vv(ii,jj))
            continue
        end
        
        % get displacements for center point
        u0 = uu(ii, jj);
        v0 = vv(ii, jj);
        
        % get linear indices of available neighbors
        rnbr = ii+roffset;
        cnbr = jj+coffset;
        keep_ind  = rnbr>=1 & rnbr<=nr & cnbr>=1 & cnbr<=nc;
        rnbr = rnbr(keep_ind); 
        cnbr = cnbr(keep_ind);        
        knbr = rnbr+(cnbr-1)*nr;
        
        % extract neighbors, keep only the num_nbr nearest non-NaN values        
        unbr = uu(knbr);
        vnbr = vv(knbr);
        %... drop NaNs
        keep_ind = ~isnan(unbr) & ~isnan(vnbr);
        unbr = unbr(keep_ind);
        vnbr = vnbr(keep_ind);
        %... keep nearest num_nbr points 
        if num_nbr < length(unbr)
            unbr = unbr(1:num_nbr);
            vnbr = vnbr(1:num_nbr);
        end
        
        % compute neighbor median, residual, and median residual 
        med_unbr = nanmedian(unbr);
        res_unbr = abs(unbr-med_unbr);
        med_res_unbr = nanmedian(res_unbr);
        
        med_vnbr = nanmedian(vnbr);
        res_vnbr = abs(vnbr-med_vnbr);
        med_res_vnbr = nanmedian(res_vnbr);
        
        % compute center normalized residual
        norm_res_u0 = abs(u0-med_unbr)/(med_res_unbr+epsilon);
        norm_res_v0 = abs(v0-med_vnbr)/(med_res_vnbr+epsilon);
        
        % combine vector components (max or sum)
        norm_res = max(norm_res_u0, norm_res_v0);
        
        % set invalid points to NaN
        if norm_res > max_norm_res
            uu(ii,jj) = NaN;
            vv(ii,jj) = NaN;
        end
        
    end
end