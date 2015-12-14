function [ini, fin, xx, yy, uu, vv] = create_dots(img_size, tform, min_spc, ...
                                          prob_white, ampl_white, ampl_black, ...
                                          sigma, max_attempts)
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
% min_spc = Scalar, minimum spacing in pixels between particles, enforced for
%   both the initial and final grids
%
% prob_white = Scalar, probability that a give particle is white
%
% ampl_white = Scalar, amplitude for gaussian model of white particles
%
% ampl_black = Scalar, amplitude for gaussian model of black particles
%
% sigma = Scalar, standard deviation of gaussian model for all particles
% 
% max_attempts = Scalar, integer, maximum number of times to attempt adding a
%   random point to the grid before considering the grid complete 
%
% %

%% initialize

% set defaults
narginchk(0,8);
if nargin == 0 || isempty(img_size)
    img_size = [100, 100]; 
end 
if nargin < 2 || isempty(tform)  
    tform = [1, 0.05,  5; 0,   1,  5]; 
end 
if nargin < 3 || isempty(min_spc)
    min_spc = 3; 
end
if nargin < 4 || isempty(prob_white)
    prob_white = 0.5;
end
if nargin < 5 || isempty(ampl_white) 
    ampl_white = 2;
end
if nargin < 6 || isempty(ampl_black)
    ampl_black = -2;
end
if nargin < 7 || isempty(sigma)
    sigma = 3;
end
if nargin < 8 || isempty(max_attempts)
    max_attempts = 1e2;
end

% check for sane inputs
validateattributes(img_size, {'numeric'}, {'integer', '>', 1, 'numel', 2});
validateattributes(tform, {'numeric'}, {'2d', 'size', [2, 3]});
validateattributes(min_spc, {'numeric'}, {'scalar', 'positive'});
validateattributes(prob_white, {'numeric'}, {'scalar', '>=' 0, '<=', 1});
validateattributes(ampl_white, {'numeric'}, {'scalar'});
validateattributes(ampl_black, {'numeric'}, {'scalar'});
validateattributes(sigma, {'numeric'}, {'scalar', 'positive'});
validateattributes(max_attempts, {'numeric'}, {'scalar', 'integer', 'positive'});
    
%% get particle locations in initial and final images
tic

% compute the reverse affine transformation of the image bounding box
x_bbox = [1, img_size(2), img_size(2),           1, 1];            
y_bbox = [1,           1, img_size(1), img_size(1), 1];              
[x_bbox_rev, y_bbox_rev] = affine_trans(tform, x_bbox, y_bbox, 0);

% get limits and footprint needed to fully populate ini and fin
xlim = [ min([x_bbox(:); x_bbox_rev(:)]); max([x_bbox(:); y_bbox_rev(:)]) ];
ylim = [ min([y_bbox(:); y_bbox_rev(:)]); max([y_bbox(:); y_bbox_rev(:)]) ];
xfoot = xlim([1, 2, 2, 1]);
yfoot = ylim([1, 1, 2, 2]);
[xfoot_fwd, yfoot_fwd] = affine_trans(tform, xfoot, yfoot, 1);

% get initial grid and its forward transform from the footprint
tri = delaunayTriangulation(xfoot, yfoot);
tri_fwd = delaunayTriangulation(xfoot_fwd, yfoot_fwd);

% add random points to grid until it is "full" (i.e. too hard to add more)
num_attempts = 0;
while num_attempts <= max_attempts
    
    % new random point and its forward transform
    xpt = rand()*range(xlim)+xlim(1);
    ypt = rand()*range(ylim)+ylim(1);
    [xpt_fwd, ypt_fwd] = affine_trans(tform, xpt, ypt, 1);
    
    % append to triangulations
    tri.Points(end+1, :) = [xpt, ypt];
    tri_fwd.Points(end+1, :) = [xpt_fwd, ypt_fwd];
    ipt = size(tri.Points, 1);
    
    % accept or reject point based on distance to neighbors 
    min_dist = min( min_dist_to_nbrs(tri, ipt), min_dist_to_nbrs(tri_fwd, ipt) );
    if min_dist >= min_spc
        % accept point
        num_attempts = 0;
    else
        % reject point
        tri.Points(ipt, :) = [];
        tri_fwd.Points(ipt, :) = [];
        num_attempts = num_attempts+1;            
    end
    
end
                    
% extract points as vectors
x_pts = tri.Points(:,1);
y_pts = tri.Points(:,2);
x_pts_fwd = tri_fwd.Points(:,1);
y_pts_fwd = tri_fwd.Points(:,2);

toc
%% generate images from particle locations
tic

% randomly sort particles into black and white colors
npts = length(x_pts);
is_white = rand(npts, 1) > prob_white ;
aa = ones(npts, 1);
aa(is_white) = ampl_white;
aa(~is_white) = ampl_black;

% init coords and images
yy = 1:img_size(1);
xx = 1:img_size(2);
ini = zeros(img_size);
fin = zeros(img_size);

% compute values at all pixel location in ini and fin
sigma2 = sigma^2;
for ii = 1:length(yy)
    for jj = 1:length(xx)
        
        vals = aa.*exp( -((x_pts-xx(jj)).^2+(y_pts-yy(ii)).^2)/sigma2 );
        ini(ii, jj) = sum(vals);
        
        vals = aa.*exp( -((x_pts_fwd-xx(jj)).^2+(y_pts_fwd-yy(ii)).^2)/sigma2 );
        fin(ii, jj) = sum(vals);
        
    end
end

toc
%% generate displacements for each pixel

[x0, y0] = meshgrid(xx, yy);
[x1, y1] = affine_trans(tform, x0(:), y0(:), 1);
x1 = reshape(x1, img_size);
y1 = reshape(y1, img_size);

uu = x1-x0;
vv = y1-y0;

%% debug

% % debug: plot grids {
% plt_xlim = [ min([x_pts; x_pts_fwd]), max([x_pts; x_pts_fwd]) ];
% plt_ylim = [ min([y_pts; y_pts_fwd]), max([y_pts; y_pts_fwd]) ];
% 
% figure
% 
% subplot(1,2,1)
% patch(x_bbox, y_bbox, 'k', 'FaceAlpha', 0.5, 'LineStyle', 'None');
% hold on
% % plot(x_pts, y_pts, 'xb'); 
% triplot(tri, 'Color', 'b');
% set(gca, 'XLim', plt_xlim, 'YLim', plt_ylim);
% 
% subplot(1,2,2)
% patch(x_bbox, y_bbox, 'k', 'FaceAlpha', 0.5, 'LineStyle', 'None');
% hold on
% % plot(x_pts_fwd, y_pts_fwd, 'xb');
% triplot(tri_fwd, 'Color', 'b');
% set(gca, 'XLim', plt_xlim, 'YLim', plt_ylim);
% % } debug

% % debug: plot images with points {
% figure
% while 1    
%     imagesc(ini);
%     hold on;
%     plot(x_pts, y_pts, 'xk');
%     title('ini');
%     hold off    
%     pause(1)
% 
%     imagesc(fin);
%     hold on;
%     plot(x_pts_fwd, y_pts_fwd, 'xk');
%     title('fin');
%     hold off
%     pause(1)
% end
% % } debug

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
nbr = unique(nbr(nbr~=idx));

% git minumum distance
dist = sqrt(sum(bsxfun(@minus, dt.Points(idx,:), dt.Points(nbr, :)).^2, 2));
min_dist = min(dist);

end