function [uu, vv] = piv_displacement(...
    ini, fin, rr, cc, uu, vv, samplen, intrlen, min_frac_data, min_frac_overlap)
% function [uu, vv] = piv_displacement(...
%     ini, fin, rr, cc, uu, vv, samplen, intrlen, min_frac_data, min_frac_overlap)
%
% Compute the displacement at initial image time on a regular grid. To
% allow for an iterative solution, the initial window locations are
% selected using "guessed" displacements, then these guesses are refined by
% the PIV calculation.
%
% Arguments:
%
%   ini, fin: 2D matrix, initial and final images at midpoint time
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
% %

% TODO: validate inputs: all grids same size, ...

% constants
min_overlap = min_frac_overlap*samplen*samplen; % frac to pixels 

% for kk = 1:numel(uu)
parfor kk = 1:numel(uu)
          
    % get sample window at initial time (no offest)
    [samp, r_samp, c_samp] = piv_window(ini, rr(kk), cc(kk), samplen);
    assert(all(size(samp) == samplen), 'Generated sample window does not match expected size');
    
    % skip if sample window is too empty
    frac_data = sum(samp(:) ~= 0)/numel(samp);
    if  frac_data < min_frac_data
        uu(kk) = NaN;
        vv(kk) = NaN;
        continue
    end
    
    % get interrogation window, offset to final time
    % NOTE: size(intr) may *not* be [intrlen, intrlen] due to rounding, it
    %   follows that the true center point may not be as requested, this is
    %   OK because we compute it directly later
    [intr, r_intr, c_intr] = piv_window(...
        fin, rr(kk) + vv(kk), cc(kk) + uu(kk), intrlen);
    
    % FIXME: strong brightess anomolies at the boundaries, these can be
    %   removed with standard AHE now
    
    % FIXME: zero-padding applied by the window picker is fucking up 
    %   the chance to use the standard normxcorr2 here.
    
    % compute masked, normalized cross correlation
    [xcr, overlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);
    
    % skip if nowhere has enough overlapping sandy pixels
    if max(overlap(:)) < min_overlap
        uu(kk) = NaN;
        vv(kk) = NaN;
        continue
    end
    
    % crop correlation plane where not enough overlapping sandy pixels
    xcr(overlap < min_overlap) = 0;
    
    % find peak with subpixel precision
    % TODO: (re)explore alternative peak-finding algorithms
    [r_peak, c_peak] = piv_peak_gauss2d(xcr);
    
    % compute displacement
    uu(kk) = c_peak - size(samp, 2) - (c_samp(1) - c_intr(1));
    vv(kk) = r_peak - size(samp, 1) - (r_samp(1) - r_intr(1));
    
    % NOTE: sample window may be only partially filled, which suggests the
    %   observation lies at the sample window centroid, not its center, we 
    %   ignore this detail here in order to keep the grid regular and avoid
    %   (more) problematic scattered interpolation
end

% report result
num_valid = sum(~isnan(uu(:)));  % same as vv
num_total = numel(uu);
fprintf('%s: valid measurements at %d/%d pts (%.2f%%)\n', ...
    mfilename, num_valid, num_total, num_valid/num_total*100);
