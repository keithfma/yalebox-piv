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
uu_fwd = xgrid_fwd-xgrid;
vv_fwd = ygrid_fwd-ygrid;

im_fwd = imwarp(im, -0.5*cat(3, uu_fwd, vv_fwd), 'cubic', 'FillValues', 0);

% find largest rectangular region that is fully populated in both im and fwd
roi = im(:,:,1)>0 & im_fwd(:,:,1)>0;

clim = [1, size(roi,2)];
for i = 1:size(roi,1)
    
    c0 = find(roi(i,:), 1, 'first');    
    if ~isempty(c0); clim(1) = max(clim(1), c0); end
    
    c1 = find(roi(i,:), 1, 'last');
    if ~isempty(c1); clim(2) = max(clim(2), c1); end
    
end

rlim = [1, size(roi,1)];
for j = 1:size(roi,2)
    
    r0 = find(roi(:,j), 1, 'first');    
    if ~isempty(r0); rlim(1) = max(rlim(1), r0); end
    
    r1 = find(roi(:,j), 1, 'last');
    if ~isempty(r1); rlim(2) = max(rlim(2), r1); end
    
end

% crop to get ini and fin
xx = xx(clim(1):clim(2));
yy = yy(rlim(1):rlim(2));
rect = [clim(1), rlim(1), clim(2)-clim(1), rlim(2)-rlim(1)];
ini = imcrop(im, rect);
fin = imcrop(im_fwd, rect);

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

