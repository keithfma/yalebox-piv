function [ini, fin, ini_roi, fin_roi, xx, yy, uu, vv] = ...
    piv_test_create_synth(img_size, tform, min_spc, prob_white, ampl_white, ...
        ampl_black, sigma, max_attempts, bnd_mean, bnd_ampl, bnd_freq)
%
% Create a synthetic image pair that consists of a random field of gaussian dots
% truncated by a boundary, which is subjected to a constant translation +
% homogenous deformation. Care is taken to generate a random grid with a known
% minimum distance between points that spans the domains of both initial and
% final images. The boundary is a sinusoidal curve, with mean, amplitude and
% frequency set by the input arguments.
%
% Input arguments:
%
% img_size = Vector, length == 2, [row, col] size of the output images
%
% tform = Matrix, size == [2, 3], affine transformation matrix using homogenous
%   coordinates. Elements include the all comoponents of the displacement
%   gradient tensor as well as constant offsets in the x and y directions:
%
%       [du/dx, du/dy, tx; 
%        dv/dx, dv/dy, ty]
%
% min_spc = Scalar, minimum spacing in pixels between particles, enforced for
%   both the initial and final (deformed) grids.
%
% prob_white = Scalar, probability that a give particlen is white. This value
% determines the "mix" of white and black sand in the synthetic images.
%
% ampl_white, ampl_black = Scalar, amplitude for gaussian model of white / black
%   particles
%
% sigma = Scalar, standard deviation of gaussian model for all particles
% 
% max_attempts = Scalar, integer, maximum number of times to attempt adding a
%   random point to the grid before considering the grid complete 
%
% bnd_mean = Scalar, mean value of the boundary sinusoid in the initial grid,
%   normalized coordinates so the range is [0, 1] show = Logical flag, 1 to show
%
% bnd_ampl = Scalar, amplitude of the boundary sinudoid in the initial grid,
%   normalized by img_size(1)
%
% bnd_freq = Scalar, frequency of the boundary sinusoid, in cycles per image
%   width
%
% Output arguments:
%
% ini, fin = 
% 
% ini_roi, fin_roi =  
% 
% xx, yy =  
% 
% uu, vv
%
% %

%% initialize

% set defaults
narginchk(0,13);
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
if nargin < 9 || isempty(bnd_mean)
    bnd_mean = 0.75;
end
if nargin < 10 || isempty(bnd_ampl)
    bnd_ampl = 0.1;
end
if nargin < 11 
    bnd_freq = 1;
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
validateattributes(bnd_mean, {'numeric'}, {'scalar'});
validateattributes(bnd_ampl, {'numeric'}, {'scalar'});
validateattributes(bnd_freq, {'numeric'}, {'scalar'});

%% get particle locations in initial and final images

% compute the reverse affine transformation of the image bounding box
x_bbox = [1, img_size(2), img_size(2),           1, 1];            
y_bbox = [1,           1, img_size(1), img_size(1), 1];              
[x_bbox_rev, y_bbox_rev] = affine_trans(tform, x_bbox, y_bbox, 0);

% get limits and footprint needed to fully populate ini and fin
xlim = [ min([x_bbox(:); x_bbox_rev(:)]); max([x_bbox(:); y_bbox_rev(:)]) ];
ylim = [ min([y_bbox(:); y_bbox_rev(:)]); max([y_bbox(:); y_bbox_rev(:)]) ];

% define boundary test function
in_bnd = @(x,y) y/img_size(1) <= ...
    bnd_mean + bnd_ampl*sin(2*pi*x/img_size(2)*bnd_freq);

% add random points to grid until it is "full" (i.e. too hard to add more)
tri = delaunayTriangulation();
tri_fwd = delaunayTriangulation();
num_attempts = 0;
while num_attempts <= max_attempts
    
    % new random point and its forward transform, skip if outside the boundary 
    xpt = rand()*range(xlim)+xlim(1);
    ypt = rand()*range(ylim)+ylim(1);    
    if ~in_bnd(xpt, ypt); continue; end    
    [xpt_fwd, ypt_fwd] = affine_trans(tform, xpt, ypt, 1);
    
    % append to triangulations
    tri.Points(end+1, :) = [xpt, ypt];
    tri_fwd.Points(end+1, :) = [xpt_fwd, ypt_fwd];
    ipt = size(tri.Points, 1);
    if ipt < 3 % build initial triangulation
        continue
    end
    
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

%% generate images from particle locations

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
        
        if in_bnd(xx(jj), yy(ii))
            vals = aa.*exp( -((x_pts-xx(jj)).^2+(y_pts-yy(ii)).^2)/sigma2 );
            ini(ii, jj) = sum(vals);
        else
            ini(ii,jj) = NaN;
        end
      
        [xx_rev, yy_rev] = affine_trans(tform, xx(jj), yy(ii), 0);
        if in_bnd(xx_rev, yy_rev)        
            vals = aa.*exp( -((x_pts_fwd-xx(jj)).^2+(y_pts_fwd-yy(ii)).^2)/sigma2 );
            fin(ii, jj) = sum(vals);
        else
            fin(ii,jj) = NaN;
        end
        
    end
end

% rescale ini and fin to the range [0, 1], and get the roi mask
min_val = min([ini(:); fin(:)]);
ini = ini-min_val;
fin = fin-min_val;

max_val = max([ini(:); fin(:)]); 
ini = ini./max_val;
fin = fin./max_val;

ini_roi = ~isnan(ini);
fin_roi = ~isnan(fin);

ini(~ini_roi) = 0;
fin(~fin_roi) = 0;

%% generate displacements for each pixel

[x0, y0] = meshgrid(xx, yy);
[x1, y1] = affine_trans(tform, x0(:), y0(:), 1);
x1 = reshape(x1, img_size);
y1 = reshape(y1, img_size);

uu = x1-x0;
vv = y1-y0;

end

%% subroutines

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