function invalid = yalebox_piv_valid_nmed(uu, vv, max_norm_res, epsilon)
% function invalid = yalebox_piv_valid_nmed(uu, vv, max_norm_res, epsilon)
%
% Validate the displacement vector field using a normalized median test with a
% 5x5 square kernel. See reference [1] for details. 
%
% Arguments:
%
%   uu, vv = 2D matrix, double, displacement vector components. NaNs are
%       treated as missing values to allow for roi masking.
%
%   max_norm_res = Scalar, double, maximum value for the normalized residual,
%       above which a vector is flagged as invalid. Reference [1] reccomends a
%       value of 2.
%
%   epsilon = Scalar, double, minumum value of the normalization factor.
%       Reference [1] recommends a value of 0.1.
%
%   invalid = 2D matrix, logical, flags identifying all invalid vectors as 1
%
%
% References: 
%
% [1] Westerweel, J., & Scarano, F. (2005). Universal outlier detection for PIV
% data. Experiments in Fluids, 39(6), 1096–1100. doi:10.1007/s00348-005-0016-6

% set defaults
if nargin < 3; max_norm_res = 2; end
if nargin < 4; epsilon = 0.1; end 

% init
[nr, nc] = size(uu);
invalid = false(nr, nc);

% % 3x3 neighborhood
% pad = 1;
% roffset = [ 1,  1,  1, ...
%             0,      0, ...
%            -1, -1, -1];
% coffset = [-1,  0,  1 ...
%            -1,      1 ...
%            -1,  0,  1];

% 5x5 neighborhood
pad = 2;
roffset = [ 2,  2,  2,  2,  2, ...
            1,  1,  1,  1,  1, ...
            0,  0,      0,  0, ...
           -1, -1, -1, -1, -1, ...
           -2, -2, -2, -2, -2];
coffset = [-2, -1,  0,  1,  2, ...
           -2, -1,  0,  1,  2, ...
           -2, -1,      1,  2, ...
           -2, -1,  0,  1,  2, ...
           -2, -1,  0,  1,  2];

% pad input arrays with NaNs to ignore out-of-bounds pixels
uu = padarray(uu, [pad, pad], NaN, 'both');
vv = padarray(vv, [pad, pad], NaN, 'both');

% loop over all positions in un-padded array
for ii = 1:nr
    for jj = 1:nc
        
        % get indices and dimensions in padded array
        iip = ii+pad;
        jjp = jj+pad;
        nrp = nr+2*pad;
        
        % get linear indices of neighbors in padded arrays
        rnbr = iip+roffset;
        cnbr = jjp+coffset;
        knbr = rnbr+(cnbr-1)*(nrp);
        
        % extract displacements for center and neighbors
        u0 = uu(iip, jjp);
        v0 = vv(iip, jjp);
        unbr = uu(knbr);
        vnbr = vv(knbr);
        
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
        
        % classify as valid or invalid
        invalid(ii, jj) = norm_res > max_norm_res;
        
    end
end

end
