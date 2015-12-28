function [ini, fin, ini_roi, fin_roi, xx, yy] = ...
    piv_test_create_image(tform, bnd_mean, bnd_ampl, bnd_freq)
%
% Create a synthetic image pair by deforming a template image with a homogenous
% deformation + constant offset, and imposing a sinusoidal boundary on the
% initial image. The image size is chosen so that both the initial and final
% (deformed) images are fully populated, and thus depends on both the
% transformation and the size of the image template.

%% initialize

% set defaults
narginchk(0,4);
if nargin == 0 || isempty(tform)  
    tform = [1, 0.05,  5; 0,   1,  5]; 
end 
if nargin < 2 || isempty(bnd_mean)
    bnd_mean = 0.75;
end
if nargin < 3 || isempty(bnd_ampl)
    bnd_ampl = 0.1;
end
if nargin < 4 || isempty(bnd_freq) 
    bnd_freq = 1;
end

% check for sane inputs
validateattributes(tform, {'numeric'}, {'2d', 'size', [2, 3]});
validateattributes(bnd_mean, {'numeric'}, {'scalar'});
validateattributes(bnd_ampl, {'numeric'}, {'scalar'});
validateattributes(bnd_freq, {'numeric'}, {'scalar'});

% local parameters
template_filename = 'test/template_fault_ss_01_sidef_251.png';

% debug: split transform into jacobian and translation
D = tform(:, 1:2);
T = tform(:, 3);

%% generate initial and final images

% load template with undeformed coordinate grid
im = imread(template_filename);
xx = 0:size(im,2)-1;
yy = 0:size(im,1)-1;
[xgrid, ygrid] = meshgrid(xx, yy);

% get forward transform displacements and image
[xgrid_fwd, ygrid_fwd] = piv_test_util_transform(tform, xgrid(:), ygrid(:), 1);
xgrid_fwd = reshape(xgrid_fwd, size(xgrid));
ygrid_fwd = reshape(ygrid_fwd, size(ygrid));
uu = xgrid_fwd-xgrid;
vv = ygrid_fwd-ygrid;

im_fwd = imwarp(im, -cat(3, uu, vv), 'cubic', 'FillValues', 0);

% find largest rectangular region that is fully populated in both im and fwd

% 
roi = im(:,:,1)>0 & im_fwd(:,:,1)>0;

[C, H, W] = FindLargestRectangles(roi, [0 0 1]);
[~, idx] = max(C(:));
[r0, c0] = ind2sub(size(roi), idx);
rect = [c0, r0, W(r0,c0), H(r0,c0)];

% crop to get ini and fin
ini = imcrop(im, rect);
fin = imcrop(im_fwd, rect);
xx = xx(rect(1):(rect(1)+rect(3)-1));
yy = yy(rect(2):(rect(2)+rect(4)-1));


%% impose boundary

% converstion from image to normalized coordinates

% generate the boundary line in normalized coordinates 
bnd_x_norm = -1:0.001:2; % 3x image width, small step
bnd_y_norm =  bnd_mean + bnd_ampl*sin(2*pi*bnd_x_norm*bnd_freq);

% convert to image coordinates
bnd_x = bnd_x_norm*range(xx)+min(xx);
bnd_y = bnd_y_norm*range(yy)+min(yy);

% get forward transform
[bnd_x_fwd, bnd_y_fwd] = piv_test_util_transform(tform, bnd_x, bnd_y, 1);

% interpolate to image coordinates and remove pixels above the boundary
bnd_y_ini = floor( interp1(bnd_x, bnd_y, xx) );
for j = 1:length(xx)
    idx = find(yy == bnd_y_ini(j));    
    ini(idx:end, j, :) = 0;
end

bnd_y_fin = floor( interp1(bnd_x_fwd, bnd_y_fwd, xx) );
for j = 1:length(xx)
    idx = find(yy == bnd_y_fin(j));    
    fin(idx:end, j, :) = 0;
end

% get roi for ini and fin
ini_roi = ini~=0;
fin_roi = fin~=0;

%% debug

% plot images
figure(1)
imagesc(rgb2gray(ini))
set(gca, 'YDir', 'normal');

figure(2)
imagesc(rgb2gray(fin))
set(gca, 'YDir', 'normal');

keyboard

function [C, H, W, M] = FindLargestRectangles(I, crit, minSize)
% finds largest rectangle regions within all points set to 1.
% input: I       - B/W boolean matrix or output of FindLargestSquares
%        minSize - [height width] - minimum width and height on regions of 
%                  interest (used to restrict final choise)
%        crit    - Optimazation Criteria parameters to optimize:
%                   crit(1)*height + crit(2)*width + crit(3)*height*width
% output: 
%         C    - value of the optimization criteria "crit" calculated for 
%                each pixel 
%         W, H - for each pixel I(r,c) return height and width of the largest 
%                all-white rectangle with its upper-left corner at I(r,c)
%         M    - Mask the largest all-white rectangle of the image
 
if (nargin<2)
  crit = [1 1 0];
end
if (nargin<3)
  minSize = [1 1];
end
p = crit;
[nR nC] = size(I);
if (minSize(1)<1), minSize(1)= floor(minSize(1)*nR); end
if (minSize(2)<1), minSize(2)= floor(minSize(2)*nC); end
if (max(I(:)) - min(I(:))==1),
  S = FindLargestSquares(I);
else
  S = I;
end
n = max(S(:));
W = S; % make a carbon copy of the matrix data
H = S;
C = ((p(1)+p(2)) + p(3)*S) .* S; % p(1)*width + p(2)*height + p(3)*height*width for height=width=S;
d = round((3*n)/4);
minH = max(min(minSize(1), d),1);
minW = max(min(minSize(2), d),1);

%% look for rectangles with width>height
hight2width = zeros(n+1,1);  % Store array with largest widths aviable for a given height
for r = 1 : nR               % each row is processed independently
  hight2width(:) = 0;        % reset the List
  for c = nC: -1 : 1         % go through all pixels in a row right to left
    s = S(r,c);              % s is a size of a square with its corner at (r,c)
    if (s>0)                 % if pixel I(r,c) is true
      MaxCrit = C(r,c);      % initialize the Max Criteria using square
      for hight = s:-1:1     % go through all possible width&hight combinations. Start with more likely to be the best
        width = hight2width(hight); % look up width for a given hight
        width = max(width+1,s);
        hight2width(hight) = width;
        Crit = p(1)*hight + p(2)*width + p(3)*width*hight;
        if (Crit>MaxCrit),   % check if it produces larger Criteria
          MaxCrit = Crit;    % if it does than save the results
          W(r,c)  = width;
          H(r,c)  = hight;
        end % if Crit
      end % for hight
      C(r,c)  = MaxCrit;
    end % if s
    hight2width((s+1):end) = 0;    % hights>s will not be aviable for the next pixel
  end % for c
end
clear hight2width

%% look for rectangles with width<height
width2hight = zeros(n+1,1);  % Store array with largest widths aviable for a given height
for c = 1 : nC               % each column is processed independently
  width2hight(:) = 0;        % reset the List
  for r = nR: -1 : 1         % go through all pixels in a column bottom to top
    s = S(r,c);              % s is a size of a square with its corner at (r,c)
    if (s>0)                 % if pixel I(r,c) is true
      MaxCrit = C(r,c);      % initialize the Max Criteria using square
      for width = s:-1:1     % go through all possible width&hight combinations. Start with more likely to be the best
        hight = width2hight(width); % look up hight for a given width
        hight = max(hight+1,s);
        width2hight(width) = hight;
        Crit = p(1)*hight + p(2)*width + p(3)*width*hight;
        if (Crit>MaxCrit),   % check if it produces larger Criteria
          MaxCrit = Crit;    % if it does than save the results
          W(r,c)  = width;
          H(r,c)  = hight;
        end % if Crit
      end % for width
      C(r,c)  = MaxCrit;
    end % if s
    width2hight((s+1):end) = 0;    % hights>s will not be aviable for the next pixel
  end % for r
end

%% Create container mask
CC = C;
CC( H<minH | W<minW ) = 0; % first try to obey size restrictions
[~, pos] = max(CC(:));
if (isempty(pos)), [~, pos] = max(C(:)); end % but when it fails than drop them
[r c] = ind2sub(size(C), pos);
M = false(size(C));
M( r:(r+H(r,c)-1), c:(c+W(r,c)-1) ) = 1;

function S = FindLargestSquares(I)
%FindLargestSquares - finds largest sqare regions with all points set to 1.
%input:  I - B/W boolean matrix
%output: S - for each pixel I(r,c) return size of the largest all-white square with its upper -left corner at I(r,c)  
[nr nc] = size(I);
S = double(I>0);
for r=(nr-1):-1:1
  for c=(nc-1):-1:1
    if (S(r,c))
      a = S(r  ,c+1);
      b = S(r+1,c  );
      d = S(r+1,c+1);
      S(r,c) = min([a b d]) + 1;
    end
  end  
end  

