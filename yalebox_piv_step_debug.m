function [] = yalebox_piv_step_debug(test_case)
% Run yalebox_piv_step for one of several debugging test cases.
%
% test_case = Scalar, integer. Select test case:
%   (1) Constant offset in x- and y-directions
%   (2) Simple shear, gamma = du/dy
%
% %

% define global parameters
create_synth = '/projectnb/glaciermod/yalebox-piv/test/synth/create';

% initialize
validateattributes(test_case, {'numeric'}, {'scalar', 'integer'});
addpath(create_synth);

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
                
    otherwise
        error('Invalid test case selected');
end

% cleanup
rmpath(create_synth);

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
