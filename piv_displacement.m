function [r1, c1, u1, v1, roi] = ...
    piv_displacement(ini, fin, r0, c0, samplen, intrlen, min_frac_data, ...
        min_frac_overlap)
%
% Compute the displacement at midpoint time from the maksed normalized cross
% correlation of sample and interrogation windows. Displacements are evaluated
% the sample grid (r0, c0), but are returned at scattered points representing
% thier position at midpoint time (r1, c1), i.e. offset by half the estimated
% displacement.
%
% Arguments:
% 
% ini, fin = 2D matrix, initial and final images at midpoint time
%
% r0, c0 = 2D matrix, row- and column-coordinates for the sample grid
%
% samplen, intrlen = Scalar, size of the sample and interrogation windows in pixels
%
% min_frac_data = Scalar, minimum fraction of the sample window that must
%   contain data (e.g. sand) for the point to be included in the ROI for PIV
%   analysis
%
% min_frac_overlap = Scalar, minimum fraction of the sample window data that
%   must overlap the interrogation window data for a point in the
%   cross-correlation to be valid
%
% r1, c1 = 2D matrix, row- and column-coordinates for the estimated displacements 
%   at midpoint time
%
% u1, v1 = 2D matrix, estimated displacements in the row (v1) and column (u1) directions 
%
% roi = 2D logical matrix, indicates whether points had enough data (e.g. sand) 
%   to estimated displacement (1) or did not (0)
% %

% constants
min_overlap = min_frac_overlap*samplen*samplen; % frac to pixels 

% allocate outputs, skipped points remain NaN, false
[nr, nc] = size(r0);
r1 = nan(nr, nc); 
c1 = nan(nr, nc); 
u1 = nan(nr, nc); 
v1 = nan(nr, nc); 
roi = true(nr, nc); 

for ii = 1:nr
    for jj = 1:nc
   
        % get sample and interrogation windows, skip if too empty
        [samp, samp_rll, samp_cll, r_centroid, c_centroid, frac_data] = ...
            piv_window(ini, r0(ii,jj), c0(ii,jj), samplen);        
        
        if frac_data < min_frac_data; 
            roi(ii,jj) = false;
            continue; 
        end
        
        [intr, intr_rll, intr_cll] = ...
            piv_window(fin, r0(ii,jj), c0(ii,jj), intrlen);
        
        % compute masked, normalized cross cor0elation
        [xcr, overlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);    
        xcr(overlap<min_overlap) = 0;
        
        % find cor0elation plane max, subpixel precision (failed pixels -> NaN)
        [rpeak, cpeak] = piv_peak_gauss2d(xcr);
        
        % convert position of the correlation max to displacement
        u1(ii,jj) = cpeak-samplen-(samp_cll-intr_cll);
        v1(ii,jj) = rpeak-samplen-(samp_rll-intr_rll);
        
        % compute location of this observation at midpoint time
        c1(ii,jj) = c_centroid + 0.5*u1(ii,jj);
        r1(ii,jj) = r_centroid + 0.5*v1(ii,jj);        
                
    end
end