function [] = yalebox_piv_step_debug(test_case)
% Run yalebox_piv_step for one of several debugging test cases.
%
% test_case = Scalar, integer. Select test case:
%   (1) Constant offset in x- and y-directions
%   (2) Simple shear, gamma = du/dy
%   (3) Wedge image pair
%
% %

% define paths to dependencies
test_synth_create = 'test/synth/create';
test_wedge_data = 'test/wedge/data';
test_wedge_prep = 'test/wedge/prep';

% initialize
validateattributes(test_case, {'numeric'}, {'scalar', 'integer'});
addpath(test_synth_create);
addpath(test_wedge_data);
addpath(test_wedge_prep);

switch test_case
    
    case 1
        
        % define case parameters
        template = 'fault_ss_01_sidef_251_template.png';
        uconst = 11; 
        vconst = 9;        
        samplen = 30;
        sampspc = 15;
        intrlen = 60;
        u0 = 0;
        v0 = 0;
        
        % create input variables
        [ini, fin, xx0, yy0, uu0, vv0] = ...
            create_constant_uv(template, uconst, vconst);
        
        % run piv
        [xx, yy, uu, vv] = ...
            yalebox_piv_step(ini, fin, xx0, yy0, samplen, sampspc, intrlen, u0, v0, 1);
        
        % analyze results
        uu_err = get_err(xx0, yy0, uu0, xx, yy, uu);
        vv_err = get_err(xx0, yy0, vv0, xx, yy, vv);
        
    case 2
        
        % define case parameters
        template = 'fault_ss_01_sidef_251_template.png';
        gamma = 0.05;
        dir = 2;
        samplen = 30;
        sampspc = 15;
        intrlen = 60;
        u0 = 0;
        v0 = 0;
        
        % create input variables
        [ini, fin, xx0, yy0, uu0, vv0] = ...
            create_simple_shear(template, gamma, dir);

        % run piv
        [xx, yy, uu, vv] = ...
            yalebox_piv_step(ini, fin, xx0, yy0, samplen, sampspc, intrlen, u0, v0, 1);
        
        % analyze results        
        uu_err = get_err(xx0, yy0, uu0, xx, yy, uu);
        vv_err = get_err(xx0, yy0, vv0, xx, yy, vv);
        
    case 3
        
        % define case parameters
        ini_file = 'fault_ss_01_sidef_250.png';
        fin_file = 'fault_ss_01_sidef_251.png';
        coord_file = 'fault_ss_01_sidef_coords.mat';
        samplen = 30;
        sampspc = 15;
        intrlen = 60;
        u0 = 0;
        v0 = 0;
        
        % prepare images
        ini = imread(ini_file);
        ini = rgb2hsv(ini);
        ini = ini(:,:,3);
        
        fin = imread(fin_file);
        fin = rgb2hsv(fin);
        fin = fin(:,:,3);
        
        % load coordinates
        load(coord_file, 'x', 'y');
        xx0 = x; clear x;
        yy0 = y; clear y;
        
        % run piv
        [xx, yy, uu, vv] = ...
            yalebox_piv_step(ini, fin, xx0, yy0, samplen, sampspc, intrlen, u0, v0, 1);
        
                
    otherwise
        error('Invalid test case selected');
end

% cleanup
rmpath(test_synth_create);
rmpath(test_wedge_data);
rmpath(test_wedge_prep);

% debug {
keyboard
% } debug

end

%% subroutines

function err = get_err(x0, y0, z0, x1, y1, z1)
%
% Interpolate exact solution to the computed solution, and take the
% difference.
%
% Arguments:
% 
%   x0, y0 = Vector, double, coordinate vectors for the exact solution
%
%   z0 = 2D matrix, double, exact solution
%
%   x1, y1 = Vector, double, coordinate vectors for the approximate
%       solution
%
%   z1 = 2D matrix, double, approximate solution
%
%   err = 2D matrix, difference between exact and approximate solutions, at
%       the resolution of the approximate solution
% %

[x1i, y1i] = meshgrid(x1, y1);
z0i = interp2(x0, y0, z0, x1i, y1i);
err = z0i-z1;

end
