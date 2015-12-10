function [ini, fin, xx, yy, uu, vv] = create_dots(img_size, A)
%
% Create a synthetic image pair that consists of a random field of gaussian dots
% truncated by a boundary, which is subjected to a constant translation +
% homogenous deformation.
%
% Arguments:
%
% img_size = Vector, length == 2, [row, col] size of the output images
%
% A = Matrix, size == [2, 3], affine transformation matrix using homogenous
%   coordinates Elements include the all comoponents of the displacement
%   gradient tensor as well as constant offsets in the x and y directions:
%       [du/dx, du/dy, tx; 
%        dv/dx, dv/dy, ty]
%
% %

%% initialize

% set defaults
narginchk(0,2);
if nargin == 0; img_size = [25, 25]; end
if nargin < 2; A = [1, 0.5, 20; 0, 1, 5]; end

% check for sane inputs
validateattributes(img_size, {'numeric'}, {'integer', '>', 1, 'numel', 2}, ...
    mfilename, 'img_size');
validateattributes(A, {'numeric'}, {'2d', 'size', [2, 3]}, mfilename, 'A');
A = [A; 0, 0, 1];

%% main

% compute the forward and reverse affine transformation of the image bounding box
bbox = [1, img_size(2), img_size(2),           1; 
        1,           1, img_size(1), img_size(1)];
bbox_homo = [bbox; ones(1,4)];
bbox_fwd_homo = A*bbox_homo;
bbox_rev_homo = inv(A)*bbox_homo;
bbox_fwd = bbox_fwd_homo(1:2, :);
bbox_rev = bbox_rev_homo(1:2, :);

% % generate a random grid
% pts = random_grid([1, img_size(2)], [1, img_size(1)], 5); % [x, y]
% pts_homo = [pts'; ones(1, size(pts, 1))]; % [x; y; ones()]
% pts_fwd_homo = A*pts_homo;
% pts_fwd = pts_fwd_homo(1:2, :)';
% pts_rev_homo = inv(A)*pts_homo;
% pts_rev = pts_rev_homo(1:2, :)';


%% debug

hold off
triplot(delaunayTriangulation(bbox'), 'Color', 'b');
hold on
triplot(delaunayTriangulation(bbox_fwd'), 'Color', 'r');
triplot(delaunayTriangulation(bbox_rev'), 'Color', 'g');
legend('ini', 'fin', 'rev')

keyboard

% dummy output arguments
ini = [];
fin = [];
xx = [];
yy = [];
uu = [];
vv = [];

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
new_pt = @() rand(1,2).*[range(xlim)+xlim(1), range(ylim)+ylim(1)];
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