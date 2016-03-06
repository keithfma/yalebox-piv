function [uu_norm, vv_norm, mm_norm, bbox] = util_normalize_displ(xx, yy, uu, vv, mm, bbox)
% 
% Normalize the input displacement vector componenets and magnitude using the
% median displacement magnitude within the specified bounding box.
%
% Magnitudes are normalised to
% facilitate comparison between steps of unequal size. This inconsistency is a
% shortcoming of the experimental apparatus. The normalization factor is a
% user-selected characteristic velocity, computed as the absolute value of the
% median of the velocity mangitude in a user-defined window. For most sandbox
% experiments, the logical choice is the incoming section.
%
% Arguments:
%   bbox = OPTIONAL
%   
%   xx, yy = 
%
%   uu, vv = 
% 
%   mm = 
%
%   uu_norm, vv_norm, mm_norm = 
% % 

if nargin < 6 || isempty(bbox)
    figure;
    imagesc(xx, yy, mm, 'AlphaData', ~isnan(mm));
    set(gca, 'YDir', 'Normal');
    title('Select bounding box for normalization using mouse');
    bbox = getrect();
    close(gcf);
end

% get normalization factor
x_in_bbox = (xx >= bbox(1)) & (xx <= (bbox(1)+bbox(3)));
y_in_bbox = (yy >= bbox(2)) & (yy <= (bbox(2)+bbox(4)));
in_bbox = logical(bsxfun(@times, x_in_bbox(:)', y_in_bbox(:)));
norm = median(mm(in_bbox));

% apply normalization
uu_norm = uu./norm;
vv_norm = vv./norm;
mm_norm = mm./norm;