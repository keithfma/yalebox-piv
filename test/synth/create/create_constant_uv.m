function [] = create_constant_uv(template, uu, vv, name)
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
%   name = String, stub for the output image pair file names, e.g. 'my_img'
%       would yield the output files my_img_1.png, my_img_2.png

% Check for sane inputs
validateattributes(template, {'char'}, {'vector'}, ...
    'create_constant_uv', 'template');
validateattributes(uu, {'numeric'}, {'scalar', 'integer'}, ...
    'create_constant_uv', 'uu');
validateattributes(uu, {'numeric'}, {'scalar', 'integer'}, ...
    'create_constant_uv', 'vv');
validateattributes(name, {'char'}, {'vector'}, ...
    'create_constant_uv', 'name');
