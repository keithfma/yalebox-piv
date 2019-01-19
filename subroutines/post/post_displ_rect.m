function [uv_median] = post_displ_rect(xx, yy, uu, vv, bbox)
%
% Compute vector quantiles of PIV displacements within a rectangular
% region. Uses vector quantile as defined in [1], rather than computing
% quantiles on each vector component separately.
%
% Arguments:
%   xx, yy: vectors, horizontal and vertical coordinates for the displacement 
%       grids uu and vv, as produced by piv_series.m 
%   uu, vv: 2D matrices, horizontal and vertical displacement grids for a
%       single time step from PIV analyses, as produced by piv_series.m
%   bbox: 4-element position vector [xmin, ymin, width, height] for the 
%       bounding box used to crop data prior to computing displacement
%       statistics
%
% Returns:
%   vector quantiles within bbox
%
% References:
% [1] Liu, Y. (2013). Noise reduction by vector median filtering. Geophysics,
%   78(3), V79–V87. http://doi.org/10.1190/geo2012-0232.1

% check input arguments
% TODO: validate remaining arguments
validateattributes(bbox, {'numeric'}, {'vector', 'numel', 4}, ...
    mfilename, 'bbox');

% extract non-NaN data within bounding box
x_in_bbox = (xx >= bbox(1)) & (xx <= (bbox(1)+bbox(3)));
y_in_bbox = (yy >= bbox(2)) & (yy <= (bbox(2)+bbox(4)));
in_bbox = logical(bsxfun(@times, x_in_bbox(:)', y_in_bbox(:)));
has_data = ~isnan(uu) & ~isnan(vv);
ubox = uu(in_bbox & has_data);
vbox = vv(in_bbox & has_data);
nbox = length(ubox);

% deal with empty (or near empty) roi
if nbox < 10;
    uv_median = [NaN, NaN];
    return
end

% compute vector statistics as defined in [1] using euclidian norm
min_norm = nan(nbox, 1);
idx = 1:nbox;
for ii = idx
   norm = sqrt((ubox(ii)-ubox(idx~=ii)).^2 + (vbox(ii)-vbox(idx~=ii)).^2);
   % norm = abs(ubox(ii)-ubox(idx~=ii) + abs(vbox(ii)-vbox(idx~=ii));
   min_norm(ii) = min(norm);
end
[~, median_idx] = min(min_norm); 
uv_median = [ubox(median_idx), vbox(median_idx)];