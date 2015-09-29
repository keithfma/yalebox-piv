function [ini, fin, xx, yy, uu, vv] = create_constant_uv(template, u, v, show)
% Create a test dataset with constant x-dir and y-dir offsets from a template
% sand image. 
%
% Case: Constant displacement in the x-direction only
%
% Arguments:
%
%   template =  String, file name of the template image, should be a subset of
%       some yalebox run, expected to be rgb
%
%   u, v = Scalar, double, constant displacements in the x and y directions,
%       in pixels, must be an integer (to avoid introducing noise by
%       interpolation).
%
%   show = Scalar, logical, flag enabling or disabling display of the test image
%       pair. If enabled, the program alternately displays the initial and final
%       images 3 times.
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
    'create_constant_uv', 'template');
validateattributes(u, {'numeric'}, {'scalar', 'integer'}, ...
    'create_constant_uv', 'u');
validateattributes(u, {'numeric'}, {'scalar', 'integer'}, ...
    'create_constant_uv', 'v');
validateattributes(show, {'numeric', 'logical'}, {'binary'}, ... 
    'create_constant_uv', 'show');

% read template data to intensity matrix
ii = imread(template);
ii = rgb2hsv(ii);
ii = ii(:, :, 3);

% apply displacements by cropping the template
if u > 0
    ini = ii(:, u+1:end);
    fin = ii(:, 1:end-u);
elseif u < 0 
    ini = ii(:, 1:end+u);
    fin = ii(:, -u+1:end);       
else
    ini = ii;
    fin = ii;    
end

if v > 0
    ini = ini(v+1:end, :);
    fin = fin(1:end-v, :);
elseif v < 0 
    ini = ini(1:end+v, :);
    fin = fin(-v+1:end, :);       
else
     % do nothing 
end

% generate pixel coordinate vectors
xx = 0:size(ini, 2)-1;
yy = 0:size(ini, 1)-1;

% generate displacement matrices
uu = u*ones(size(ini));
vv = v*ones(size(ini));

% display image pair
if show
    clr = [min(ii(:)), max(ii(:))];
    h = figure;
    for i = 1:3
        figure(h)
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






