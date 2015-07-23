function [imgc] = correct_apply(img, mask, baseline, scale)
%
% Apply mask and color corrections computed using the correct_mask() and
% correct_color functions.
%
% Arguments:
%
%   img =
%   mask = 
%   baseline = 
%   scale =
%   imgc = 
%
% Keith Ma, July 2015

%% check for sane inputs

%% apply corrections

imgc = img;
imgc(repmat(~mask, [1, 1, 3])) = 0;
imgc = bsxfun(@minus, imgc, baseline);
imgc = bsxfun(@times, imgc, scale);
imgc = min(1, max(0, imgc));
