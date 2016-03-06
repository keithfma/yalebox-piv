function [uu_norm, vv_norm, mm_norm] = ...
    plot_normalize_displacement(bbox, xx, yy, uu, vv, mm)
% 
% Normalize the input displacement vector componenets and magnitude using the
% median displacement magnitude within the specified bounding box.
%
% Arguments:
%   bbox = 
%   
%   xx, yy = 
%
%   uu, vv = 
% 
%   mm = 
%
%   uu_norm, vv_norm, mm_norm = 
% % 

% get normalization factor
x_in_bbox = (xx >= bbox(1)) & (xx <= (bbox(1)+bbox(3)));
y_in_bbox = (yy >= bbox(2)) & (yy <= (bbox(2)+bbox(4)));
in_bbox = logical(bsxfun(@times, x_in_bbox(:)', y_in_bbox(:)));
norm = median(mm(in_bbox));

% apply normalization
uu_norm = uu./norm;
vv_norm = vv./norm;
mm_norm = mm./norm;