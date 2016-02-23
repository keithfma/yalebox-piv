function [rr, cc, uu, vv, roi] = piv_displacement(xcr, r0, c0, u0, v0)
%
% Compute the displacement from the cross correlation plane.

% preallocate outputs, skipped points retain these values ([] or NaN)
[nr, nc] = size(r0);
rr = nan(nr, nc); % point location at midpoint time
cc = nan(nr, nc); % point location at midpoint time
uu = nan(nr, nc);
vv = nan(nr, nc);
roi = true(nr, nc);

% find displacements from correlation planes
for ii = 1:nr
    for jj = 1:nc
        
        % skip if no xcr is available
        if isempty(xcr{ii,jj})
            roi(ii,jj) = false;
            continue;
        end
   
        % find correlation plane max, subpixel precision (failed pixels -> NaN)
        [rpeak, cpeak] = piv_peak_gauss2d(xcr{ii,jj});
        
        % convert position of the correlation max to displacement
        uu(ii,jj) = cpeak + u0(ii,jj);
        vv(ii,jj) = rpeak + v0(ii,jj);
        
        % compute location of this observation at midpoint time
        cc(ii,jj) = c0(ii,jj) + 0.5*uu(ii,jj);
        rr(ii,jj) = r0(ii,jj) + 0.5*vv(ii,jj);        
        
    end
end
