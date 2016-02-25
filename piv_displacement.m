function [r_pts, c_pts, u_pts, v_pts, roi] = ...
    piv_displacement(ini, fin, rr, cc, samplen, intrlen)
% function [r_pts, c_pts, u_pts, v_pts, roi] = ...
%     piv_displacement(ini, fin, rr, cc, samplen, intrlen)
%
% Compute the displacement at midpoint time from the maksed normalizd cross
% correlation of sample and interrogation windows.
%
% Arguments:
% 
% ini, fin = 
% rr, cc = 
% samplen, intrlen = 
% r_pts, c_pts = 
% u_pts, v_pts = 
% roi = 
% %

% constants
min_frac_data = 0.5;
min_frac_overlap = min_frac_data/2;
min_overlap = min_frac_overlap*samplen*samplen; % frac to pixels 

% allocate outputs, skipped points remain NaN, false
[nr, nc] = size(rr);
r_pts = nan(nr, nc); 
c_pts = nan(nr, nc); 
u_pts = nan(nr, nc); 
v_pts = nan(nr, nc); 
roi = true(nr, nc); 

for ii = 1:nr
    for jj = 1:nc
   
        % get sample and interrogation windows, skip if too empty
        [samp, samp_rll, samp_cll, r_centroid, c_centroid, frac_data] = ...
            piv_window(ini, rr(ii,jj), cc(ii,jj), samplen);        
        
        if frac_data < min_frac_data; 
            roi(ii,jj) = false;
            continue; 
        end
        
        [intr, intr_rll, intr_cll] = ...
            piv_window(fin, rr(ii,jj), cc(ii,jj), intrlen);
        
        % compute masked, normalized cross correlation
        [xcr, overlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);    
        xcr(overlap<min_overlap) = 0;
        
        % find correlation plane max, subpixel precision (failed pixels -> NaN)
        [rpeak, cpeak] = piv_peak_gauss2d(xcr);
        
        % convert position of the correlation max to displacement
        u_pts(ii,jj) = cpeak-samplen-(samp_cll-intr_cll);
        v_pts(ii,jj) = rpeak-samplen-(samp_rll-intr_rll);
        
        % compute location of this observation at midpoint time
        c_pts(ii,jj) = c_centroid + 0.5*u_pts(ii,jj);
        r_pts(ii,jj) = r_centroid + 0.5*v_pts(ii,jj);        
                
    end
end