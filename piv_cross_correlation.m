function [xcr, r0, c0, u0, v0] = ...
             piv_cross_correlation(ini, fin, rr, cc, samplen, intrlen, ...
                min_frac_data,  min_frac_overlap)
%
% Compute masked normalized cross correlation for all sample points. Returns the
% correlation planes, the centroids of the sample windows, and the offset needed
% to convert peak position to dispacement.

%   u_offset, v_offset = Scalar, double, constant offset between the location of
%       the correlation plane peak and displacement, such that u = c_peak +
%       u_offset, and v = r_peak + v_offset. Accounts for padding in
%       cross-correlation (-len) and also for relative position of interogation 
%       and sample windows (e.g.for columns: -(samp_min_col - intr_min_col))

min_overlap = min_frac_overlap*samplen*samplen; % frac to pixels 

% preallocate outputs, skipped points retain these values ([] or NaN)
[nr, nc] = size(rr);
xcr = cell(nr, nc);
r0 = nan(nr, nc); % sample window centroid row
c0 = nan(nr, nc); % sample window centroid column
u0 = nan(nr, nc); % displacement offset, u = c_peak + u0
v0 = nan(nr, nc); % displacement offset, v = r_peak + v0

for ii = 1:nr
    for jj = 1:nc
   
        [samp, samp_rll, samp_cll, r0(ii,jj), c0(ii,jj), frac_data] = ...
            piv_window(ini, rr(ii,jj), cc(ii,jj), samplen);
        
        if frac_data < min_frac_data; continue; end
        
        [intr, intr_rll, intr_cll] = ...
            piv_window(fin, rr(ii,jj), cc(ii,jj), intrlen);
        
        u0(ii,jj) = -samplen-(samp_cll-intr_cll);
        v0(ii,jj) = -samplen-(samp_rll-intr_rll);
        
        [xcr{ii,jj}, overlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);    
        xcr{ii,jj}(overlap<min_overlap) = 0;
        
    end
end