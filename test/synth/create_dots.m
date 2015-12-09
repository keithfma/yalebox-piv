function [ini, fin, xx, yy, uu, vv] = create_dots()
%
% Create a synthetic image pair that consists of a random field of gaussian dots
% truncated by a boundary, which is subjected to a constant translation +
% homogenous deformation.
%
% Arguments:
%
% %

% debug: dummy input arguments {
% } debug 

pts = random_grid(1000, [0, 1], [0, 1], 0.01);

% debug: dummy output arguments {
ini = [];
fin = [];
xx = [];
yy = [];
uu = [];
vv = [];
% } debug

function pts = random_grid(num_pts, xlim, ylim, min_spc)
%
% Generate a set of num_pts random points within the limits xlim and ylim with a
% minimum spacing of min_spc.
%
% Arguments:
%
% num_pts = Scalar, integer, number of points in the output set
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
num_attempts = 0;

% add random points until there are num_pts of them in the set
while size(tri.Points, 1) < num_pts && num_attempts <= max_num_attempts
    
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
