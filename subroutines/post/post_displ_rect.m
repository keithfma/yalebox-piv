function [uv_median] = post_displ_rect(xx, yy, uu, vv, bbox, show)
%
% Compute vector quantiles of PIV displacements within a rectangular
% region. Uses vector quantile as defined in [1], rather than computing
% quantiles on each vector component separately.
%
% Arguments:
% 
%   xx, yy: vectors, horizontal and vertical coordinates for the displacement 
%       grids uu and vv, as produced by piv_series.m 
% 
%   uu, vv: 2D matrices, horizontal and vertical displacement grids for a
%       single time step from PIV analyses, as produced by piv_series.m
% 
%   bbox: 4-element position vector [xmin, ymin, width, height] for the 
%       bounding box used to crop data prior to computing displacement
%       statistics
% 
%   show: set true to plot the the bounding box and related information, 
%       default is false
%
% Returns:
%   vector quantiles within bbox
%
% References:
% [1] Liu, Y. (2013). Noise reduction by vector median filtering. Geophysics,
%   78(3), V79�V87. http://doi.org/10.1190/geo2012-0232.1

% set defaults
narginchk(5, 6);
if nargin == 5; show = false; end

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

if nbox >= 10
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
    
else
    % deal with empty (or near empty) roi
    uv_median = [NaN, NaN];
    
end

% display is requested
if show
    figure
    
    ax1 = subplot(2, 2, 1:2);
    im1 = imagesc([xx(1), xx(end)], [yy(1), yy(end)], sqrt(uu.^2 + vv.^2));
    im1.AlphaData = ~isnan(uu);
    hold on;
    plot([bbox(1), bbox(1) + bbox(3), bbox(1) + bbox(3), bbox(1)          , bbox(1)], ...
         [bbox(2), bbox(2)          , bbox(2) + bbox(4), bbox(2) + bbox(4), bbox(2)], ...
         '-k', 'LineWidth', 1);
    legend('bbox');
    axis equal tight
    ax1.Title.String = 'Displacement Magnitude';
    ax1.XAxis.Label.String = 'x [m]';
    ax1.YAxis.Label.String = 'y [m]';
    ax1.YDir = 'normal';
    
    ax2 = subplot(2, 2, 3);
    histogram(ax2, ubox);
    hold on;
    plot([uv_median(1), uv_median(1)], ax2.YLim, 'r', 'LineWidth', 2);
    legend('distribution', 'vector median')
    ax2.Title.String = 'U in Bounding Box';
    ax2.XAxis.Label.String = 'U [m/step]';
    ax2.YAxis.Label.String = 'Count';
    
    ax3 = subplot(2, 2, 4);
    histogram(ax3, vbox);
    hold on;
    plot([uv_median(2), uv_median(2)], ax3.YLim, 'r', 'LineWidth', 2);
    legend('distribution', 'vector median')
    ax3.Title.String = 'V in Bounding Box';
    ax3.XAxis.Label.String = 'V [m/step]';
    ax3.YAxis.Label.String = 'Count';
end