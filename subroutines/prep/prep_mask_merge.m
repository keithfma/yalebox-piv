function [mask, mask_bnd] = prep_mask_merge(mask_auto, mask_manual, show)
% function [mask, mask_bnd] = prep_mask_merge(mask_auto, mask_manual, show)
%
% Merge automatic and manual masks, and return the merged mask and its
%   boundary. Note that the boundary is smoothed to provide a more stable
%   estimate.
%
% Arguments:
%   mask_auto: 2D matrix, logical, the automatic mask returned by
%       mask_prep_auto
%   mask_manual: 2d matrix, logical, the manually-defined mask returned by
%       mask_prep_manual
%   show: set true to display results, or false to run quietly
%
% Returns:
%   mask: 2D matrix, logical, the combined mask after boundary smoothing
%   bnd.top: [row, col] positions of the top of the mask
%   bnd.bottom: [row, col] positions of the bottom of the mask
%   bnd.left: [row, col] positions of the left side of the mask
%   bnd.right: [row, col] positions of the right side of the mask
%   bnd.complete: [row, col] positions of the complete mask (in-order
%       merge of the top, bottom, left, and right components
% %

narginchk(2, 3);
if nargin < 3; show = false; end

validateattributes(mask_auto, {'logical'}, {'2d'});
validateattributes(mask_manual, {'logical'}, {'2d', 'size', size(mask_auto)});

% combine the masks before smoothing boundaries
mask_combined = mask_auto & mask_manual;

% FIXME: will need to deal with the "pits" generated by singe masking
% upstream, the smoothing looks about right, but cannot account for these,
% perhaps a morphological filter? or a heuristic on finding outliers?

% FIXME: STOP - saving this boundary adds extra work, let's just leave
%   it as a mask, and use this code only when we need it.

% extract the boundaries and smooth them
% note: there are tricky edge cases here that prevent a simpler approach,
%   smoothing the whole boundaries as paramteric curves will work for side
%   and top views, single and multipart sand bodies, etc.
raw_mask_bnds = bwboundaries(mask_combined);
mask_bnds = cell(size(raw_mask_bnds));
for ii = 1:length(raw_mask_bnds)
    npts = 3;
    rc = raw_mask_bnds{ii};
    r = rc(:, 1);
    c = rc(:, 2);
    delta = sqrt(diff(r).^2 + diff(c).^2);
    dist = [0; cumsum(delta)];
    sr = smooth(dist, r, npts/length(r), 'lowess');
    sc = smooth(dist, c, npts/length(c), 'lowess');
    mask_bnds{ii} = [sr, sc];    
end

% create a new pixel mask from the smoothed boundaries

mask = [];
mask_bnd = [];

if show
    disp('display a plot');
end