
function [r_tm, c_tm, u_tm, v_tm, roi] = ...
    piv_displacement_shift(ini, fin, r_tm, c_tm, u_tm, v_tm, samplen, intrlen, min_frac_data, ...
        min_frac_overlap, high_quality, verbose)
%
% NEW VERSION: displace sample and interrogation windows from midpoint time,
% assumes NO image deformation. Below docs are out of date
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
% verbose = Scalar, display verbose messages (1) or don't (0)
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
[nr, nc] = size(r_tm);
roi = false(nr, nc); 

for ii = 1:nr
    for jj = 1:nc
   
        % get sample window, offset to initial time
        % NOTE: size(samp) may *not* be [samplen, samplen] due to rounding
        % NOTE: it follows that the center point must be re-computed
        r_ti = r_tm(ii, jj) - 0.5*v_tm(ii, jj);
        c_ti = c_tm(ii, jj) - 0.5*u_tm(ii, jj);
        [samp, r_samp, c_samp] = piv_window_shift(ini, r_ti, c_ti, samplen);        
        r_ti = r_samp(1) + 0.5*(r_samp(end) - r_samp(1));
        c_ti = c_samp(1) + 0.5*(c_samp(end) - c_samp(1));
        
        % skip if sample window is too empty
        frac_data = sum(samp(:) ~= 0)/numel(samp);
        if  frac_data < min_frac_data; 
            continue; 
        end        
        
        % get interrogation window, offset to final time
        % NOTE: size(intr) may *not* be [intrlen, intrlen] due to rounding
        % NOTE: it follows that the center point must be re-computed
        r_tf = r_tm(ii, jj) + 0.5*v_tm(ii, jj);
        c_tf = c_tm(ii, jj) + 0.5*u_tm(ii, jj);
        [intr, r_intr, c_intr] = piv_window_shift(fin, r_tf, c_tf, intrlen);
        r_tf = r_intr(1) + 0.5*(r_intr(end) - r_intr(1));
        c_tf = c_intr(1) + 0.5*(c_intr(end) - c_intr(1));
        
        % compute masked, normalized cross correlation
        [xcr, overlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);
        
        % skip if nowhere has enough overlapping sandy pixels
        if max(overlap(:)) < min_overlap
            continue
        end
        
        % crop correlation plane where not enough overlapping sandy pixels
        xcr(overlap < min_overlap) = 0;
        
        % find peak with subpixel precision
        
        % TODO: include two options, a 'rough' option for initial passes, and a
        % 'fine' option for the last pass
        
        % NOTE: two options below have approx the same accuracy in aggregate,
        %   the first is ~25% faster, but error distribution is multimodal,
        %   second is slower, but error distribution is approx normal.
        if high_quality
            [r_peak, c_peak] = piv_peak_optim_interp(xcr, 1e-6); % FINE
        else
            [r_peak, c_peak] = piv_peak_gauss2d(xcr); % ROUGH
        end
        
        % compute relative displacement, subpixel precision
        dc = c_peak - size(samp, 2) - (c_samp(1) - c_intr(1));
        dr = r_peak - size(samp, 1) - (r_samp(1) - r_intr(1));
        
        % compute total displacement
        u_tm(ii, jj) = c_tf - c_ti + dc;
        v_tm(ii, jj) = r_tf - r_ti + dr;

        % get centroid of the sample window
        [r_idx, c_idx] = find(samp ~= 0);
        num = length(r_idx);
        r0 = sum(r_samp(r_idx))/num;
        c0 = sum(c_samp(c_idx))/num;
        
        % update observation point to midpoint time
        r_tm(ii, jj) = r0 + 0.5*v_tm(ii, jj);
        c_tm(ii, jj) = c0 + 0.5*u_tm(ii, jj);
        
        % add to ROI
        roi(ii, jj) = true;

        % % DEBUG: line added to pause routine and generate plots, if desired
        % if all(intr(:)~=0)
        %     save('piv_dump.mat');
        %     error('SAVED WINDOW DATA TO piv_dump.mat');
        % end
                
    end
end

% set points outside the ROI to NaN
u_tm(~roi) = NaN;
v_tm(~roi) = NaN;

if verbose
    num_total = numel(roi);
    num_found = sum(~isnan(u_tm(:)));
    num_skipped = sum(~roi(:));
    num_dropped = sum(roi(:) & isnan(u_tm(:)));
    fprintf('%s: found %d (%.1f%%), skipped %d (%.1f%%), dropped %d (%.1f%%)\n', ...
        mfilename, num_found, num_found/num_total*100, num_skipped, ...
        num_skipped/num_total*100, num_dropped, num_dropped/num_total*100);
end