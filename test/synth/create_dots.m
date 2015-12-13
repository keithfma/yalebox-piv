function [ini, fin, xx, yy, uu, vv] = create_dots(img_size, tform)
%
% Create a synthetic image pair that consists of a random field of gaussian dots
% truncated by a boundary, which is subjected to a constant translation +
% homogenous deformation.
%
% Arguments:
%
% img_size = Vector, length == 2, [row, col] size of the output images
%
% tform = Matrix, size == [2, 3], affine transformation matrix using homogenous
%   coordinates Elements include the all comoponents of the displacement
%   gradient tensor as well as constant offsets in the x and y directions:
%       [du/dx, du/dy, tx; 
%        dv/dx, dv/dy, ty]
%
% %

%% initialize

% debug parameters
min_spc = 2;

% set defaults
narginchk(0,2);
if nargin == 0; 
    img_size = [25, 25]; 
end
if nargin < 2;  
    tform = [1, 0.5, 20; 
             0,   1,  5]; 
end

% check for sane inputs
validateattributes(img_size, {'numeric'}, {'integer', '>', 1, 'numel', 2}, ...
    mfilename, 'img_size');
validateattributes(tform, {'numeric'}, {'2d', 'size', [2, 3]}, mfilename, 'tform');

%% main

% compute the reverse affine transformation of the image bounding box
x_bbox = [1, img_size(2), img_size(2),           1, 1];            
y_bbox = [1,           1, img_size(1), img_size(1), 1];              
[x_bbox_rev, y_bbox_rev] = affine_trans(tform, x_bbox, y_bbox, 0);

% get the footprint of the points needed to fully populate ini and fin
pts_xlim = [ min([x_bbox(:); x_bbox_rev(:)]); max([x_bbox(:); y_bbox_rev(:)]) ];
pts_ylim = [ min([y_bbox(:); y_bbox_rev(:)]); max([y_bbox(:); y_bbox_rev(:)]) ];

% generate an intial dense grid and its forward transform
[xx, yy] = meshgrid(pts_xlim(1):pts_xlim(2), pts_ylim(1):pts_ylim(2));
xx = xx(:); 
yy = yy(:);
[xx_fwd, yy_fwd] = affine_trans(tform, xx, yy, 1);

% triangulate the initial and deformed grids (for finding neighbors)
tri = delaunayTriangulation([xx, yy]);
tri_fwd = delaunayTriangulation([xx_fwd, yy_fwd]);

% point vectors are not used again
clear xx yy xx_fwd yy_fwd
    
% prune initial grid to desired spacing
npts = length(tri.Points(:,1));
order = randperm(npts);
for ii = 1:npts % iterate over all points in random order
    
    % get minimum spacing in original and deformed grids
    kk = order(ii);
    d0 = min_dist_to_nbrs(tri, kk);
    d1 = min_dist_to_nbrs(tri_fwd, kk);
    
    % remove point if it fall below the minimum spacing
    if min(d0, d1) < min_spc
        tri.Points(kk, :) = [];
        tri_fwd.Points(kk ,:) = [];
        order(order>kk) = order(order>kk)-1;
    end
    
end

% extract points and delete triangulations
x_pts = tri.Points(:,1);
y_pts = tri.Points(:,2);
clear tri
x_pts_fwd = tri_fwd.Points(:,1);
y_pts_fwd = tri_fwd.Points(:,2);
clear tri_fwd

%% debug

% debug: plot grids {
plt_xlim = [ min([x_pts; x_pts_fwd]), max([x_pts; x_pts_fwd]) ];
plt_ylim = [ min([y_pts; y_pts_fwd]), max([y_pts; y_pts_fwd]) ];

figure

subplot(1,2,1)
patch(x_bbox, y_bbox, 'k', 'FaceAlpha', 0.5, 'LineStyle', 'None');
hold on
plot(x_pts, y_pts, 'xb'); 
set(gca, 'XLim', plt_xlim, 'YLim', plt_ylim);

subplot(1,2,2)
patch(x_bbox, y_bbox, 'k', 'FaceAlpha', 0.5, 'LineStyle', 'None');
hold on
plot(x_pts_fwd, y_pts_fwd, 'xb'); 
set(gca, 'XLim', plt_xlim, 'YLim', plt_ylim);
% } debug

keyboard

% dummy output arguments
ini = [];
fin = [];
xx = [];
yy = [];
uu = [];
vv = [];

end

function [x_pts_out, y_pts_out] = affine_trans(tform, x_pts_in, y_pts_in, fwd)
% function [x_pts_out, y_pts_out] = affine_trans(tform, x_pts_in, y_pts_in, fwd)
% 
% tform = 2x3 affine transformation matrix
%
% x_pts_in, y_pts_in = vectors of x, y point coordinates
%
% fwd = Scalar, logical, flag indicating if the transform should be forward (1)
%   or reverse (0)
%
% x_pts_out, y_pts_out = vectors of transformed x,y point coordinates
%
% %

% get full transform matrix
A = [tform; 0 0 1]; 
if ~fwd; 
    A = inv(A); 
end

% transform points, maintaining vector shape
pts_in = [x_pts_in(:)'; y_pts_in(:)'; ones(1, length(x_pts_in))];

pts_out = A*pts_in;

x_pts_out = reshape(pts_out(1,:), size(x_pts_in));
y_pts_out = reshape(pts_out(2,:), size(y_pts_in));

end

function min_dist = min_dist_to_nbrs(dt, idx)
% function min_dist = min_dist_to_nbrs(dt, idx)
%
% find minimum distance to neighboring vertices in the delaunay triangulation
% 'dt' for the vertex with index 'idx'
%
% %

 % get indices of all connected neighbors from triangulation 
attach = cell2mat(vertexAttachments(dt, idx));
nbr = dt.ConnectivityList(attach, :);
nbr = unique(nbr(:));
nbr = nbr(nbr~=idx);

% git minumum distance
dist = sqrt(sum(bsxfun(@minus, dt.Points(idx,:), dt.Points(nbr, :)).^2, 2));
min_dist = min(dist);

end


% function pts = random_grid_old(xlim, ylim, min_spc)
% %
% % Generate a set of num_pts random points within the limits xlim and ylim with a
% % minimum spacing of min_spc.
% %
% % Arguments:
% %
% % xlim, ylim = Vectors, length==2 , [minimum, maximum] coordinates for points in
% %   the set
% %
% % min_spc = Scalar, minimum permissible distance between any pair of points in
% %   the set
% %
% % %
% 
% % parameters
% max_num_attempts = 1e3;
% 
% % initialize
% new_pt = @() rand(1,2).*[range(xlim), range(ylim)]+[xlim(1), ylim(1)];
% tri = delaunayTriangulation([new_pt(); new_pt(); new_pt()]);
% 
% % add random points until there are num_pts of them in the set
% num_attempts = 0;
% while num_attempts <= max_num_attempts
%     
%     % insert a new point into triangulation
%     tri.Points(end+1, :) = new_pt();
%     
%     % find distance to all connected neighbors for the new point
%     attach = cell2mat(vertexAttachments(tri, size(tri.Points, 1)));
%     nbr = tri.ConnectivityList(attach, :);
%     nbr = unique(nbr(:));
%     nbr = nbr(nbr~=size(tri.Points, 1));
%     dist = sqrt(sum(bsxfun(@minus, tri.Points(end,:), tri.Points(nbr, :)).^2, 2));
%     
%     % accept or reject the new point
%     if min(dist) >= min_spc
%         num_attempts = 0;
%     else
%         tri.Points(end, :) = [];
%         num_attempts = num_attempts+1;            
%     end
%     
% end
% 
% pts = tri.Points;
% 
% end