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

%% generate initial and final images

% load template and convert to normalized grayscale 
im = imread(template_filename);
im = rgb2hsv(im);
im = im(:,:,3);
im = prep_intensity(im, true(size(im)), 31);

% get undeformed coordinate grid
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

% crop to largest rectangular region that is fully populated in both im and fwd
%... uses 3rd party functions from the file exchange, see private folder

roi = im(:,:,1)>0 & im_fwd(:,:,1)>0;

if any(roi == 0)
    
    [C, H, W] = FindLargestRectangles(roi, [0 0 1]);
    [~, idx] = max(C(:));
    [r0, c0] = ind2sub(size(roi), idx);
    rect = [c0, r0, W(r0,c0), H(r0,c0)];
    
    % crop to get ini and fin
    ini = imcrop(im, rect);
    fin = imcrop(im_fwd, rect);
    xx = xx(rect(1):(rect(1)+rect(3)-1));
    yy = yy(rect(2):(rect(2)+rect(4)-1));
    
else
    
    ini = im;
    fin = im_fwd;
end



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

% %% debug
% 
% % plot images
% figure(1)
% imagesc(rgb2gray(ini))
% set(gca, 'YDir', 'normal');
% 
% figure(2)
% imagesc(rgb2gray(fin))
% set(gca, 'YDir', 'normal');

