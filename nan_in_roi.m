function [] = nan_in_roi(data, roi)
% 
% Check for any NaNs in the ROI of data, enter debugger if found.
% %

if any(isnan(data(roi)))
    fprintf('%s: Found NaN in ROI\n', mfilename);
    keyboard
end
