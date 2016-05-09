function img_a = movie_aux_spoint(img, x, y, tip, len, clr, opac)
% function img_a = movie_aux_spoint(img, x, y, tip, len, clr, opac)
%
% Superimpose s-point triangles on the input image. Wrapper function for
% the MATLAB function insertShape from the Computer Vision System Toolbox.
%
% Arguments:
%
%   img_a = 3D matrix, the output annotated image
%
%   img = 2D or 3D matrix, the input image to be annotated 
%
%   x, y = Vector, double,  coordinate vectors for the input image "img"
%
%   tip = N x 2 matrix, double, where N is the number of triangles to draw.
%       Each row gives the location of the triangle tip as x, y in world
%       coordinates
%
%   len = Scalar, double, equilateral triangle side length in world coordinates
%
%   clr = MATLAB color definintion, can be a valid string or RGB triplet
%
%   opac = Scalar, double, normalized opacity in the range [0, 1]
%
% Keith Ma, August 2015

% check for sane inputs

% initialize
nx = numel(x);
ny = numel(y);
nt = size(tip, 1); % 1 triangle per row
dx = abs(x(1)-x(2));
dy = abs(y(1)-y(2));
pos = nan(nt, 6); 

% define position matrix
for i = 1:nt
   rr = interp1(y, 1:ny, tip(i,2), 'linear', 'extrap'); 
   cc = interp1(x, 1:nx, tip(i,1), 'linear', 'extrap');
      
   if abs(rr-ny) < abs(rr) % tip at top
       pos(i,:) = [cc,          rr,                  ...
                   cc+len/2/dx, rr-sqrt(3)/2*len/dy, ...
                   cc-len/2/dx, rr-sqrt(3)/2*len/dy];       
   
   else % tip at bottom
       pos(i,:) = [cc,          rr,                  ...
                   cc+len/2/dx, rr+sqrt(3)/2*len/dy, ...
                   cc-len/2/dx, rr+sqrt(3)/2*len/dy];
   end   
end

% annotate image    
img_a = insertShape(img, 'FilledPolygon', pos, ...
        'Color', clr, ...
        'Opacity', opac);
    
% done