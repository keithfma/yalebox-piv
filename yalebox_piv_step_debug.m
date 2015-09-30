function [] = yalebox_piv_step_debug(test_case)
% Run yalebox_piv_step for one of several debugging test cases.
%
% test_case = Scalar, integer. Select test case:
%   (1) Constant offset in x- and y-directions
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
        uconst = 10; 
        vconst = 10;        
        samplen = 30;
        sampspc = 15;
        intrlen = 60;
        u0 = 0;
        v0 = 0;
        
        % create input variables
        [ini, fin, xx0, yy0, uu0, vv0] = ...
            create_constant_uv('fault_ss_01_sidef_251_template.png', uconst, vconst);
        
        % run piv
        [xx, yy, uu, vv] = ...
            yalebox_piv_step(ini, fin, xx0, yy0, samplen, sampspc, intrlen, u0, v0, 1);
        
        % analyze results
                
    otherwise
        error('Invalid test case selected');
end

% cleanup
rmpath(create_synth);