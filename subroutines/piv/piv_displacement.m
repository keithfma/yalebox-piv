function [r_tm, c_tm, u_tm, v_tm, roi] = piv_displacement(...
    ini, fin, r_tm_i, c_tm_i, u_tm_i, v_tm_i, samplen, intrlen, ...
    min_frac_data, min_frac_overlap, high_quality)
%
% Compute the displacement at midpoint time from the maksed normalized cross
% correlation of sample and interrogation windows. Displacements are evaluated
% at the specified points (r_tm, c_tm), but are returned at (scattered) points
% representing thier position at midpoint time (an updated r_tm, c_tm), i.e.
% offset by half the estimated displacement. To allow for an iterative solution,
% the initial window locations are selected using "guessed" displacements, then
% these guesses are refined by the PIV calculation.
%
% Arguments:
%
% ini, fin: 2D matrix, initial and final images at midpoint time
% r_tm_i, c_tm_i: 2D matrix, row- and column-coordinates for the initial sample grid
% u_tm_i, v_tm_i: 2D matrix, initial guess for x (column) and y (row)
%   displacements at points in r_tm_i, c_tm_i.
% samplen, intrlen = Scalar, size of the sample and interrogation windows in pixels
% min_frac_data = Scalar, minimum fraction of the sample window that must
%   contain data (e.g. sand) for the point to be included in the ROI for PIV
%   analysis
% min_frac_overlap = Scalar, minimum fraction of the sample window data that
%   must overlap the interrogation window data for a point in the
%   cross-correlation to be valid
%
% r_tm, c_tm = Vector, row- and column-coordinates for the estimated
%   displacements at midpoint time
% u_tm, v_tm = Vector, estimated displacements in the x (column) and y (row)
%   directions
% %

% constants
min_overlap = min_frac_overlap*samplen*samplen; % frac to pixels 

% allocate outputs, skipped points remain NaN, false
[nr, nc] = size(r_tm_i);
nn = numel(r_tm_i);
u_tm = nan(nr, nc);
v_tm = nan(nr, nc);
r_tm = nan(nr, nc);
c_tm = nan(nr, nc);

parfor kk = 1:nn
          
    % get sample window, offset to initial time
    % NOTE: size(samp) may *not* be [samplen, samplen] due to rounding, it
    %   follows that the true center point may not be r_ti, c_ti, that's OK
    r_ti = r_tm_i(kk) - 0.5*v_tm_i(kk);
    c_ti = c_tm_i(kk) - 0.5*u_tm_i(kk);
    [samp, r_samp, c_samp] = piv_window(ini, r_ti, c_ti, samplen);
    
    % skip if sample window is too empty
    frac_data = sum(samp(:) ~= 0)/numel(samp);
    if  frac_data < min_frac_data;
        continue;
    end
    
    % get interrogation window, offset to final time
    % NOTE: size(intr) may *not* be [intrlen, intrlen] due to rounding, it
    %   follows that the true center point may not be r_ti, c_ti, that's OK
    r_tf = r_tm_i(kk) + 0.5*v_tm_i(kk);
    c_tf = c_tm_i(kk) + 0.5*u_tm_i(kk);
    [intr, r_intr, c_intr] = piv_window(fin, r_tf, c_tf, intrlen);
    
    % compute masked, normalized cross correlation
    [xcr, overlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);
    
    % skip if nowhere has enough overlapping sandy pixels
    if max(overlap(:)) < min_overlap
        continue
    end
    
    % crop correlation plane where not enough overlapping sandy pixels
    xcr(overlap < min_overlap) = 0;
    
    % find peak with subpixel precision
    % NOTE: two options below have approx the same accuracy in aggregate,
    %   the "low-quality" option is faster, but has a multimodal error
    %   distribution, the "high-quality" option is slower, but yields a
    %   approx normal error distribution.
    if high_quality
        [r_peak, c_peak] = piv_peak_optim_interp(xcr, 1e-6);
    else
        [r_peak, c_peak] = piv_peak_gauss2d(xcr);
    end
    
    % compute displacement
    u_tm(kk) = c_peak - size(samp, 2) - (c_samp(1) - c_intr(1));
    v_tm(kk) = r_peak - size(samp, 1) - (r_samp(1) - r_intr(1));
    
    % get centroid of the sample window
    [r_idx, c_idx] = find(samp ~= 0);
    num = length(r_idx);
    r_samp_cntr = sum(r_samp(r_idx))/num;
    c_samp_cntr = sum(c_samp(c_idx))/num;
    
    % update observation point to midpoint time
    r_tm(kk) = r_samp_cntr + 0.5*v_tm(kk);
    c_tm(kk) = c_samp_cntr + 0.5*u_tm(kk);
end

% convert matrices to vectors of valid measurements only
roi = ~isnan(u_tm); % same as v_tm, r_tm, c_tm
u_tm = u_tm(roi);
v_tm = v_tm(roi);
r_tm = r_tm(roi);
c_tm = c_tm(roi);

% report result
num_valid = numel(u_tm);
num_total = numel(roi);
fprintf('%s: valid measurements at %d/%d pts (%.2f%%)\n', ...
    mfilename, num_valid, num_total, num_valid/num_total*100);
