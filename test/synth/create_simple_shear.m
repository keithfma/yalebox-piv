function [ini, fin, xx, yy, uu, vv] = create_simple_shear(template, gamma, dir, show)
% Create a test dataset with simple shear deformation from a template sand
% image.
%
% Arguments:
%
%   template =  String, file name of the template image, should be a subset of
%       some yalebox run, expected to be rgb
%
%   gamma = shear component of the displacement gradient matrix, if dir is 1,
%       then gamma is du/dy, if dir is 2, then gamma is dv/dx, units are
%       pixels/pixel
%
%   dir = Scalar, integer, flag indicating the direction to apply simple shear,
%       1 for the y-direction, and 2 for the x direction.
%
%   show = Scalar, logical, flag enabling or disabling display of the test image
%       pair. If enabled, the program alternately displays the initial and final
%       images 5 times.
%
%   ini, fin = 2D matrix, double, initial and final intensity images
%
%   xx, yy = Vector, double, coordinate vectors for ini and fin, in pixels
%
%   uu, vv = 2D matrix, double, exact displacements for this test, in pixels
% %
 
% set defaults
if nargin < 4
    show = false;
end

% check for sane inputs
validateattributes(template, {'char'}, {'vector'}, ...
    'create_simple_shear', 'template');
validateattributes(gamma, {'numeric'}, {'scalar'}, ...
    'create_simple_shear', 'gamma');
assert(gamma ~= 0, 'gamma must be non-zero');
validateattributes(dir, {'numeric'}, {'scalar', 'integer'}, ...
    'create_simple_shear', 'dir');
assert(ismember(dir, [1, 2]), 'dir must be either 1 or 2');
validateattributes(show, {'numeric', 'logical'}, {'binary'}, ... 
    'create_simple_shear', 'show');

% read template data to normalized intensity matrix
ini = imread(template);
ini = rgb2hsv(ini);
ini = yalebox_prep_intensity(ini, true(size(ini(:,:,1))), 31);

% compute displacement fields
[nr, nc] = size(ini);
[xx, yy] = meshgrid(0:nc-1, 0:nr-1);
if dir == 1
    % dv/dx = gamma
    uu = zeros(nr, nc);
    vv = gamma*xx;
elseif dir == 2
    % du/dy = gamma
    uu = gamma*yy;
    vv = zeros(nr, nc);
end

% apply transformation, trim range overshoot due to interpolation
fin = imwarp(ini, cat(3, -uu, -vv), 'cubic', 'FillValues', 0);
fin(fin<0) = 0+eps;
fin(fin>1) = 1-eps;

% trim images and displacements
if dir == 1

    % drop cols that have all nans
    ckeep = ~all(isnan(fin), 1);
    fin = fin(:, ckeep);
    ini = ini(:, ckeep);  
    uu = uu(:, ckeep);    
    vv = vv(:, ckeep);  
    
    % drop rows that have any nans
    rkeep = ~any(isnan(fin), 2);
    fin = fin(rkeep, :);
    ini = ini(rkeep, :);
    uu = uu(rkeep, :);
    vv = vv(rkeep, :);
        
elseif dir == 2
    
    % drop rows that have all nans
    rkeep = ~all(isnan(fin), 2);
    fin = fin(rkeep, :);
    ini = ini(rkeep, :);
    uu = uu(rkeep, :);
    vv = vv(rkeep, :);
    
    % drop cols that have any nans
    ckeep = ~any(isnan(fin), 1);
    fin = fin(:, ckeep);
    ini = ini(:, ckeep);  
    uu = uu(:, ckeep);    
    vv = vv(:, ckeep);    
   
end

% generate pixel coordinate vectors
xx = 0:size(ini, 2)-1;
yy = 0:size(ini, 1)-1;

% display image pair
if show
    clr = [min(ini(:)), max(ini(:))];
    h = figure;
    for i = 1:5
        figure(h);
        imagesc(ini);
        set(gca, 'YDir', 'normal')        
        set(gcf, 'Name', 'ini');
        caxis(clr);
        pause(1)
        imagesc(fin); 
        set(gca, 'YDir', 'normal')
        set(gcf, 'Name', 'fin');
        caxis(clr);
        pause(1)
    end
end
