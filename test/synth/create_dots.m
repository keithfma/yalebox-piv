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

bbox = [          1,           1;
                  1, img_size(2);
        img_size(1), img_size(2);
        img_size(1),           1;
                  1,           1];
              
bbox_rev = affine_transform(tform, bbox, 0);

% get the footprint of the points needed to fully populate ini and fin
pts_xlim = [ min([bbox(:,1); bbox_rev(:,1)]); max([bbox(:,1); bbox_rev(:,1)]) ];
pts_ylim = [ min([bbox(:,2); bbox_rev(:,2)]); max([bbox(:,2); bbox_rev(:,2)]) ];

% generate a random grid
pts = random_grid(pts_xlim, pts_ylim, 5); % [x, y]


%% debug

% plot point grids 
hold off
plot(bbox(:,1), bbox(:,2), '-k');
hold on
plot(bbox_rev(:,1), bbox_rev(:,2), '--k');
% triplot(delaunayTriangulation(pts), 'Color', 'r');
plot(pts(:,1), pts(:,2), 'xr');
legend('img', 'rev', 'pts')

keyboard

% dummy output arguments
ini = [];
fin = [];
xx = [];
yy = [];
uu = [];
vv = [];

end

function pts_out = affine_transform(tform, pts_in, fwd)
% function pts_out = affine_transform(tform, pts_in, fwd)
% 
% tform = 2x3 affine transformation matrix
%
% pts_in = Nx2 matrix of [x, y] point coordinates
%
% fwd = Scalar, logical, flag indicating if the transform should be forward (1)
%   or reverse (0)
%
% pts_out = Nx2 matrix of transformed [x,y] point coordinates
% %

% get full transform matrix
A = [tform; 0 0 1]; 
if ~fwd; 
    A = inv(A); 
end

% transform
pts_in = [pts_in, ones(size(pts_in, 1), 1)]';
pts_out = A*pts_in;
pts_out = pts_out(1:2,:)';

end

function pts = random_grid(xlim, ylim, min_spc)
%
% Generate a set of num_pts random points within the limits xlim and ylim with a
% minimum spacing of min_spc.
%
% Arguments:
%
% xlim, ylim = Vectors, length==2 , [minimum, maximum] coordinates for points in
%   the set
%
% min_spc = Scalar, minimum permissible distance between any pair of points in
%   the set
%
% %

% parameters
max_num_attempts = 1e3;

% initialize
new_pt = @() rand(1,2).*[range(xlim), range(ylim)]+[xlim(1), ylim(1)];
tri = delaunayTriangulation([new_pt(); new_pt(); new_pt()]);

% add random points until there are num_pts of them in the set
num_attempts = 0;
while num_attempts <= max_num_attempts
    
    % insert a new point into triangulation
    tri.Points(end+1, :) = new_pt();
    
    % find distance to all connected neighbors for the new point
    attach = cell2mat(vertexAttachments(tri, size(tri.Points, 1)));
    nbr = tri.ConnectivityList(attach, :);
    nbr = unique(nbr(:));
    nbr = nbr(nbr~=size(tri.Points, 1));
    dist = sqrt(sum(bsxfun(@minus, tri.Points(end,:), tri.Points(nbr, :)).^2, 2));
    
    % accept or reject the new point
    if min(dist) >= min_spc
        num_attempts = 0;
    else
        tri.Points(end, :) = [];
        num_attempts = num_attempts+1;            
    end
    
end

pts = tri.Points;

end