function [translate_row, translate_col, quality] = piv_register(...
    sample_image, sample_mask, sample_origin, ...
    reference_image, reference_mask, reference_origin, ...
    min_overlap ...
)
%
% Compute the [row, column] translation needed to register the sample to the reference image
%
% Arguments:
%   sample_image: single-band image (2D array) to be registered
%   sample_mask: boolean array, true where sample_image has valid data, and false elsewhere
%   sample_origin: the intrinsic coordinates, [row, column], of sample_image(1, 1)
%   reference_image: single-band image (2D array) to register sample against
%   reference_mask: boolean array, true where reference_image has valid data, and false elsewhere
%   reference_origin: the intrinsic coordinates, [row, column], of reference_image(1, 1)
%   min_overlap: minimum number of pixels of the sample image data that must overlap the reference
%       image data for a point in the cross-correlation to be valid
%
% Returns:
%   delta_row, delta_col: scalar floats, transalation to apply to get register the sample_image
%       to reference_image as best we can
%   quality: enumerated quality flag used to communicate any known quality issues / failures 
%       encountered during registration
%
% Feature flags:
%   YALEBOX_PIV_REGISTER_METHOD: select which registration method to use. This flag is for 
%       development / debugging only, and will be removed when the methodology is settled. 
%       Available options are:
%           + 'normxcorr2_masked' [default]
%           + more to come soon!
% %

% feature flags
method = feature_flag('YALEBOX_PIV_REGISTER_METHOD', 'normxcorr2_masked');


translate_row = NaN;  % will return null translation if registration fails
translate_col = NaN; 


if strcmp(method, 'normxcorr2_masked')
    % Masked Normalized Cross-Correlation
    
    % compute masked, normalized cross correlation
    [xcr, overlap] = normxcorr2_masked(reference_image, sample_image, reference_mask, sample_mask);
        
    % fail if nowhere has enough overlapping sandy pixels
    if max(overlap(:)) < min_overlap
        quality = Quality.BelowMinOverlap;
        return
    end
    
    % crop correlation plane where not enough overlapping sandy pixels
    xcr(overlap < min_overlap) = 0;
    
    % find peak with subpixel precision
    % TODO: (re)explore alternative peak-finding algorithms
    [r_peak, c_peak] = piv_peak_gauss2d(xcr);
    if isnan(r_peak) || isnan(c_peak)
        quality = Quality.PeakFindingFailed;
        return
    end
 
    % compute displacement
    translate_row = r_peak - size(sample_image, 1) - (sample_origin(1) - reference_origin(1));
    translate_col = c_peak - size(sample_image, 2) - (sample_origin(2) - reference_origin(2));
    quality = Quality.Valid;
    
    
else
    % Bad value for method
    error('Unknown method specified: %s', method);

end