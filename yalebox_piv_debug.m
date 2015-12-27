function [xx, yy, uu, vv] = yalebox_piv_debug(test_case)
% Run yalebox_piv for one of several debugging test cases.
%
% test_case = Scalar, integer. Select test case:
%   (1) Constant offset in x- and y-directions
%   (2) Simple shear, gamma = du/dy
%   (3) Wedge image pair
%   (4) Synthetic image, simple shear with offset
%
% %

% setup environment
addpath('test/synth');
addpath('test/wedge');

switch test_case
    
    case 1
        
        % case parameters
        template = 'template_fault_ss_01_sidef_251.png';
        uconst = 11; 
        vconst = 9;        
        
        % piv parameters
        samplen = [60, 30];
        sampspc = [30, 15];
        intrlen = [120, 60];
        npass = [3, 3];
        valid_max = 2;
        valid_eps = 0.01;
        
        % create input variables
        [ini, fin, xx0, yy0, uu0, vv0] = ...
            create_constant_uv(template, uconst, vconst);
        
        % run piv
        [xx, yy, uu, vv] = yalebox_piv(ini, fin, ini~=0, fin~=0, xx0, yy0, samplen, ...
                               sampspc, intrlen, npass, valid_max, ...
                               valid_eps, 1);
        
        % analyze results
        uu_err = get_err(xx0, yy0, uu0, xx, yy, uu);
        vv_err = get_err(xx0, yy0, vv0, xx, yy, vv);
        print_err_qnt(uu_err, vv_err);
        show_err(uu_err, vv_err);
        
    case 2
        
        % case parameters
        template = 'template_fault_ss_01_sidef_251.png';
        gamma = 0.05;
        dir = 1;
        
        % piv parameters
        samplen = [30];%, 20];
        sampspc = [15];%, 10];
        intrlen = [60];%, 30];
        npass = [3];%, 5];      
        valid_max = 2;
        valid_eps = 0.01;
        
        % create input variables
        [ini, fin, xx0, yy0, uu0, vv0] = ...
            create_simple_shear(template, gamma, dir);

        % run piv
        [xx, yy, uu, vv] = yalebox_piv(ini, fin, ini~=0, fin~=0, xx0, yy0, samplen, ...
                               sampspc, intrlen, npass, valid_max, ...
                               valid_eps, 1);
                           
        % analyze results        
        uu_err = get_err(xx0, yy0, uu0, xx, yy, uu);
        vv_err = get_err(xx0, yy0, vv0, xx, yy, vv);
        print_err_qnt(uu_err, vv_err);
        show_err(uu_err, vv_err);
        
    case 3
        
        % case parameters
        data_file = 'fault_ss_01_sidef_250_251.mat';
%         data_file = 'fault_ss_01_sidef_025_026.mat';
        
        % piv parameters
        samplen = [30, 30];%, 15];
        sampspc = [30, 15];%, 15];
        intrlen = [120, 60];%, 30];
        npass = [1, 3];%, 3];
        valid_max = 2;
        valid_eps = 0.01;
        
        % load wedge data
        load(data_file, 'xx', 'yy', 'ini', 'fin');
        xx0 = xx; %#ok!
        yy0 = yy; %#ok!
        
        % run piv
        [xx, yy, uu, vv] = yalebox_piv(ini, fin, ini~=0, fin~=0, xx0, yy0, samplen, ...
                               sampspc, intrlen, npass, valid_max, ...
                               valid_eps, 1);  
                    
        % analyze results
        [xxg, yyg] = meshgrid(xx, yy);                           
        [displacement,spin,Dv,Dd,D2x,D2y,WkStar,AkStar] = ...
            yalebox_decompose_step(xxg, yyg, uu, vv, ~isnan(uu));
        
        figure('units', 'normalized', 'position', [0.05, 0.05, 0.9, 0.9]);
        
        subplot(2,1,1)
        imagesc(xx, yy, displacement);
        colorbar
        hold on
        quiver(xx, yy, uu, vv, 2, '-k');
        axis equal
        axis tight        
        hold off
        
        subplot(2,1,2)
        imagesc(xx, yy, Dd)
        colorbar
        axis equal
        axis tight        
        
    case 4
        
        % case parameters
        img_size = [500, 500];
        tform = [1, 0.05,  0; ...
                 0,    1,  0];
        min_spc = 3;
        prob_white = 0.5;
        ampl_white = 2;
        ampl_black = -2;
        sigma = 3;
        max_attempts = 1e2;
        bnd_mean = 1.1;
        bnd_ampl = 0;
        bnd_freq = 1;
        
        init_from_file = 1;
        init_filename = 'test/synth/init_synth.mat';
        
        % piv parameters
        samplen = [30];%, 30];
        sampspc = [30];%, 15];
        intrlen = [100];%, 60];
        npass = [1];%, 3];      
        valid_max = 2;
        valid_eps = 0.01;
        
        % init input variables, either from file (fast) or computationally (slow)
        if init_from_file && exist(init_filename, 'file')==2 
            vars = load(init_filename, 'ini', 'fin', 'ini_roi', 'fin_roi', ...
                'xx0', 'yy0', 'uu0', 'vv0');
            ini = vars.ini;
            fin = vars.fin;
            ini_roi = vars.ini_roi;
            fin_roi= vars.fin_roi;
            xx0 = vars.xx0;
            yy0 = vars.yy0;
            uu0 = vars.uu0;
            vv0 = vars.vv0;
            
        else
            [ini, fin, ini_roi, fin_roi, xx0, yy0, uu0, vv0] = ...
                create_dots(img_size, tform, min_spc, prob_white, ampl_white, ...
                    ampl_black, sigma, max_attempts, bnd_mean, bnd_ampl, ...
                    bnd_freq, 1);
                
            save(init_filename, 'ini', 'fin', 'ini_roi', 'fin_roi', 'xx0', ...
                    'yy0', 'uu0', 'vv0');
        end
                   
        % run piv
        [xx, yy, uu, vv] = ...
            yalebox_piv(ini, fin, ini_roi, fin_roi, xx0, yy0, samplen, ...
                sampspc, intrlen, npass, valid_max, valid_eps, 1);
                           
        % analyze results        
        uu_err = get_err(xx0, yy0, uu0, xx, yy, uu);
        vv_err = get_err(xx0, yy0, vv0, xx, yy, vv);
        print_err_qnt(uu_err, vv_err);
        show_err(uu_err, vv_err);
                
    otherwise
        error('Invalid test case selected');
end

% cleanup environment
addpath('test/synth');
addpath('test/wedge');

% % debug {
% keyboard
% % } debug

end

%% subroutines

function err = get_err(x0, y0, z0, x1, y1, z1)
% function err = get_err(x0, y0, z0, x1, y1, z1)
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

function [] = print_err_qnt(uerr, verr)
%
% Print a table of absolute error quantiles 
%
% Arguments:
% 
%   uerr, verr = Double, error matrix, difference between exact and approximate
%       solutions for displacement in the x- and y-directions
% %   

qnt = 0:0.10:1;
qu = quantile(abs(uerr(:)), qnt);
qv = quantile(abs(verr(:)), qnt);
qm = quantile(sqrt(uerr(:).^2+verr(:).^2), qnt);

fprintf('displacement vector error quantiles\n');
fprintf('qnt\t| uu\t\t| vv\t\t| mag\n'); 
for i = 1:length(qnt)
    fprintf('%.2f\t| %.2e\t| %.2e\t| %.2e\n', ...
        qnt(i), qu(i), qv(i), qm(i));
end
fprintf('\n');

end

function [] = show_err(uerr, verr)
%
% Display a quick plot of the error matrices
 
figure('units', 'normalized', 'position', [0.05, 0.3, 0.9, 0.5]);

subplot(1, 3, 1)
imagesc(uerr);
set(gca, 'YDir', 'normal');
colorbar
title('uu error')

subplot(1, 3, 2)
imagesc(verr);
set(gca, 'YDir', 'normal');
colorbar
title('vv error')

subplot(1, 3, 3)
imagesc(sqrt(uerr.^2+verr.^2));
set(gca, 'YDir', 'normal');
colorbar
title('error magnitude')

end

