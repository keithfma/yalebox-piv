function [] = piv_test_run_wedge()
%
% Run a hard-coded test case using a images from a wedge experiment.
%
% %

%% parameters

% piv parameters
samplen = [30, 30];%, 15];
sampspc = [30, 15];%, 15];
intrlen = [120, 60];%, 30];
npass = [1, 3];%, 3];
valid_max = 2;
valid_eps = 0.01;

% image parameters
dir_name = '~/Documents/dissertation/yalebox-exp-fault/data/fault_ss_01/image/crop_sidef/';
ini_file = 'fault_ss_01_sidef_250.png';
fin_file = 'fault_ss_01_sidef_251.png';
hue_lim = [0, 0.5]; 
val_lim = [0, 0.5];
entr_lim = [0.5, 1];
entr_win = 9;
morph_open_rad = 15;
morph_erode_rad = 10;
eql_nwin = 31;

% local parameters
data_file = 'test/wedge.mat';

%% generate (or load) test images

% check if defined parameters match saved parameters
try
    F = load(data_file, 'dir_name', 'ini_file', 'fin_file', 'hue_lim', ...
        'val_lim', 'entr_lim', 'entr_win', 'morph_open_rad', ...
        'morph_erode_rad', 'eql_nwin');
    same = ...
        strcmp(F.dir_name, dir_name) && ...
        strcmp(F.ini_file, ini_file) && ...
        strcmp(F.fin_file, fin_file) && ...
        all(F.hue_lim == hue_lim) && ...
        all(F.val_lim == val_lim) && ...
        all(F.entr_lim == entr_lim) && ...
        F.entr_win == entr_win && ...
        F.morph_open_rad == morph_open_rad && ...
        F.morph_erode_rad == morph_erode_rad && ...
        F.eql_nwin == eql_nwin;
catch
    same = 0;
end

% load PIV input data
if same
    fprintf('parameters are not modified, loading input variables from file\n');
    F = load(data_file, 'ini', 'ini_roi', 'fin', 'fin_roi', 'xx', 'yy');
    ini = F.ini;
    ini_roi = F.ini_roi;
    fin = F.fin;
    fin_roi = F.fin_roi;
    xx = F.xx;
    yy = F.yy;
    
else
    fprintf('Parameters are modified, generating new input variables\n');
    [ini, fin, ini_roi, fin_roi, xx, yy] = ...
        piv_test_create_wedge(dir_name, ini_file, fin_file, hue_lim, val_lim, ...
            entr_lim, entr_win, morph_open_rad, morph_erode_rad, eql_nwin);
    save(data_file, 'dir_name', 'ini_file', 'fin_file', 'hue_lim', ...
        'val_lim', 'entr_lim', 'entr_win', 'morph_open_rad', ...
        'morph_erode_rad', 'eql_nwin', 'ini', 'ini_roi', 'fin', 'fin_roi', ...
        'xx', 'yy');
end

clear F

%% run PIV and analyze results

% run piv
[xx, yy, uu, vv] = yalebox_piv(ini, fin, ini_roi, fin_roi, xx, yy, samplen, ...
    sampspc, intrlen, npass, valid_max, valid_eps, 1);

% compute strain values
[xgrid, ygrid] = meshgrid(xx, yy);
[displ, ~, ~, Dd, ~, ~, ~, ~] = ...
    yalebox_decompose_step(xgrid, ygrid, uu, vv, ~isnan(uu));

% print and plot standard results
piv_test_util_plot(xx, yy, uu, vv, displ, Dd);