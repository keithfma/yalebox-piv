function [ini, fin, xx, yy] = create_constant_uv(template, uu, vv, show)
% Create a test dataset with constant x-dir and y-dir offsets from a template
% sand image. Outputs are cropped, grayscale versions of the template image.
%
% Case: Constant displacement in the x-direction only
%
% Arguments:
%
%   template_file =  String, file name of the template image, should be a subset
%       of some yalebox run, expected to be rgb
%
%   uu, vv = Scalar, double, constant displacements in the x and y directions,
%       in pixels, must be an integer (to avoid introducing noise by
%       interpolation).
%
%   show = Scalar, logical, flag enabling or disabling display of the test image
%       pair. If enabled, the program enters an infinite loop, alternately
%       displaying one image and the other until the user terminates via the
%       Ctrl+C command or closes the figure window.
% %

% check for sane inputs
validateattributes(template, {'char'}, {'vector'}, ...
    'create_constant_uv', 'template');
validateattributes(uu, {'numeric'}, {'scalar', 'integer'}, ...
    'create_constant_uv', 'uu');
validateattributes(uu, {'numeric'}, {'scalar', 'integer'}, ...
    'create_constant_uv', 'vv');
validateattributes(show, {'numeric', 'logical'}, {'binary'}, ... 
    'create_constant_uv', 'show');

% read template data to intensity matrix
ii = imread(template);
ii = rgb2hsv(ii);
ii = ii(:, :, 3);

% apply displacements by cropping the template
if uu > 0
    ini = ii(:, uu+1:end);
    fin = ii(:, 1:end-uu);
elseif uu < 0 
    ini = ii(:, 1:end+uu);
    fin = ii(:, -uu+1:end);       
else
    ini = ii;
    fin = ii;    
end

if vv > 0
    ini = ini(vv+1:end, :);
    fin = fin(1:end-vv, :);
elseif vv < 0 
    ini = ini(1:end+vv, :);
    fin = fin(-vv+1:end, :);       
else
     % do nothing 
end

% display image pair
if show
    clr = [min(ii(:)), max(ii(:))];
    figure;
    while 1
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






