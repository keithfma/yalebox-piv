function [uu, vv, qual] = piv_displacement(...
    ini, fin, rr, cc, uu, vv, samplen, intrlen, min_frac_data, min_frac_overlap)
% function [uu, vv, qual] = piv_displacement(...
%     ini, ini_mask, fin, fin_mask, rr, cc, uu, vv, samplen, intrlen, min_frac_data, min_frac_overlap)
%
% Compute the displacement at initial image time on a regular grid. To
% allow for an iterative solution, the initial window locations are
% selected using "guessed" displacements, then these guesses are refined by
% the PIV calculation.
%
% Arguments:
%
%   ini, fin: 2D matrix, initial and final images, expected to be floats with NaN anywhere there
%       is no data
% 
%   rr, cc: 2D matrix, row- and column-coordinates for the initial sample grid
%
%   uu, vv: 2D matrix, initial guess for x (column) and y (row)
%       displacements at points in rr, cc.
%
%   samplen, intrlen = Scalar, size of the sample and interrogation windows in pixels
%
%   min_frac_data = Scalar, minimum fraction of the sample window that must
%       contain data (e.g. sand) for the point to be included in the ROI for PIV
%       analysis
%
%   min_frac_overlap = Scalar, minimum fraction of the sample window data that
%       must overlap the interrogation window data for a point in the
%       cross-correlation to be valid
%
% Returns:
%
%   uu, vv = Matrix, estimated displacements in the x (column) and y (row)
%       directions
%   qual: Matrix, Quality flags indicating how each observation was
%       computed (or not computed), possible values are enumerated in the
%       Quality class
% %

% TODO: validate inputs: all grids same size, ...

% FIXME: with mirror pad, may be no need to test for full
%   sample window or overlap fraction, but this requires first fixing
%   the zero-padding applied by sample windows. Once this is done, we can
%   also remove the quality flag from this function

% FIXME: not sure I am handling edges right...  

% constants
min_overlap = min_frac_overlap*samplen*samplen; % frac to pixels

% switch to using 0 to signify missing data (to avoid numerical oblivion),
%   and store logical mask to safely recover this info (since 0 is a valid pixel value)
ini_mask = ~isnan(ini);  % true where has data
ini(~ini_mask) = 0;

fin_mask = ~isnan(fin);
fin(~fin_mask) = 0;

% allocate quality and mask matrices
qual = repmat(Quality.Valid, size(uu));

for kk = 1:numel(uu)
    
    % get sample window at initial time (no offest)
    [samp, r_samp, c_samp] = piv_window(ini, rr(kk), cc(kk), samplen);
    [samp_mask, ~, ~] = piv_window(ini_mask, rr(kk), cc(kk), samplen);
    assert(all(abs(size(samp) - samplen) <= 1), ...  # allows 1-pixel leeway, favor centered window over exact specified size
        'Generated sample window does not match expected size');
    
    % skip if sample window is too empty
    % FIXME: is this needed when we use mirror padding?
    frac_data = sum(samp_mask(:))/numel(samp);
    if  frac_data < min_frac_data
        qual(kk) = Quality.EmptySampleWindow;
        continue
    end
    
    % get interrogation window, offset to final time
    % NOTE: size(intr) may *not* be [intrlen, intrlen] due to rounding, it
    %   follows that the true center point may not be as requested, this is
    %   OK because we compute it directly later
    [intr, r_intr, c_intr] = piv_window(...
        fin, rr(kk) + vv(kk), cc(kk) + uu(kk), intrlen);
    [intr_mask, ~, ~] = piv_window(...
        fin_mask, rr(kk) + vv(kk), cc(kk) + uu(kk), intrlen);
    
    % compute masked, normalized cross correlation
    try
        [xcr, overlap] = normxcorr2_masked(intr, samp, intr_mask, samp_mask);
    catch
        keyboard
    end
        
    % skip if nowhere has enough overlapping sandy pixels
    if max(overlap(:)) < min_overlap
        qual(kk) = Quality.BelowMinOverlap;
        continue
    end
    
    % crop correlation plane where not enough overlapping sandy pixels
    xcr(overlap < min_overlap) = 0;
    
    % find peak with subpixel precision
    % TODO: (re)explore alternative peak-finding algorithms
    [r_peak, c_peak] = piv_peak_gauss2d(xcr);
    if isnan(r_peak) || isnan(c_peak)
        qual(kk) = Quality.PeakFindingFailed;
        continue
    end
    
    % compute displacement
    uu(kk) = c_peak - size(samp, 2) - (c_samp(1) - c_intr(1));
    vv(kk) = r_peak - size(samp, 1) - (r_samp(1) - r_intr(1));
    
    % NOTE: sample window may be only partially filled, which suggests the
    %   observation lies at the sample window centroid, not its center, we 
    %   ignore this detail here in order to keep the grid regular and avoid
    %   (more) problematic scattered interpolation
end

% report result
num_valid = sum(qual(:) == Quality.Valid);
num_total = numel(uu);
fprintf('%s: valid measurements at %d/%d pts (%.2f%%)\n', ...
    mfilename, num_valid, num_total, num_valid/num_total*100);
