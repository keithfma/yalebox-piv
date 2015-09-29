function [ini, fin, xx, yy] = create_simple_shear(template, gamma, dir, show)
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
%       pair. If enabled, the program enters an infinite loop, alternately
%       displaying one image and the other until the user terminates via the
%       Ctrl+C command or closes the figure window.
%
%   ini, fin = 2D matrix, double, initial and final intensity images
%
%   xx, yy = Vector, double, coordinate vectors for ini and fin, in pixels
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

% read template data to intensity matrix
ini = imread(template);
ini = rgb2hsv(ini);
ini = ini(:, :, 3);

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

% apply transformation
fin = imwarp(ini, cat(3, -uu, -vv), 'cubic', 'FillValues', NaN);

% trim both initial and final images
if dir == 1
    if gamma > 0
        i0 = find(~isnan(fin(:, end)), 1, 'first');
        i1 = size(fin, 1);
    elseif gamma <0
        i0 = 1;
        i1 = find(~isnan(fin(:, end)), 1, 'last');
    end
    ini = ini(i0:i1, :);
    fin = fin(i0:i1, :);
        
elseif dir == 2
    if gamma > 0 
        j0 = find(~isnan(fin(end, :)), 1, 'first');
        j1 = size(fin, 2);
    elseif gamma < 0
        j0 = 1;
        j1 = find(~isnan(fin(end, :)), 1, 'last');
        % do something
    end
    ini = ini(:, j0:j1);
    fin = fin(:, j0:j1);
end

% generate pixel coordinate vectors
xx = 0:size(ini, 2)-1;
yy = 0:size(ini, 1)-1;

% display image pair
if show
    clr = [min(ini(:)), max(ini(:))];
    h = figure;
    while 1
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
