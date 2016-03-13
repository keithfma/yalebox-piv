function [] = test_piv_run_wedge(force)
%
% Run a hard-coded test case using a images from a wedge experiment.
%
% Arguments:
%
% force = Scalar, binary, flag indicating whether to force recomputing input
%   variables (1) or not (0). This is useful if the downstream code has changed,
%   but input parameters have not.
% %

%% parameters

% piv parameters
samplen = [30, 30];
sampspc = [15, 15];
intrlen = [120, 60];
npass = [1, 4];
valid_max = 2;
valid_eps = 0.1;
lowess_span_pts = 16;
spline_tension = 0.95;
min_frac_data = 0.8;
min_frac_overlap = min_frac_data/2;

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
func_name = 'test_piv_run_wedge';

%% parse arguments and set defaults

narginchk(0,1);
if nargin == 0 || isempty(force) 
    force = 0;
end

validateattributes(force, {'numeric'}, {'scalar', 'binary'});

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

% report status
if same 
    fprintf('%s: Parameters are not modified\n', func_name); 
else
    fprintf('%s: Parameters are modified\n', func_name); 
end
if force
    fprintf('%s: Force recompute enabled\n', func_name);
else
    fprintf('%s: Force recompute disabled\n', func_name);
end

% load PIV input data
if same && ~force
    fprintf('%s: Loading input variables from file\n', func_name);
    F = load(data_file, 'ini', 'ini_roi', 'fin', 'fin_roi', 'xx', 'yy');
    ini = F.ini;
    ini_roi = F.ini_roi;
    fin = F.fin;
    fin_roi = F.fin_roi;
    xx = F.xx;
    yy = F.yy;
    
else
    fprintf('%s: Generating new input variables\n', func_name);
    [ini, fin, ini_roi, fin_roi, xx, yy] = ...
        test_piv_create_wedge(dir_name, ini_file, fin_file, hue_lim, val_lim, ...
            entr_lim, entr_win, morph_open_rad, morph_erode_rad, eql_nwin);
    save(data_file, 'dir_name', 'ini_file', 'fin_file', 'hue_lim', ...
        'val_lim', 'entr_lim', 'entr_win', 'morph_open_rad', ...
        'morph_erode_rad', 'eql_nwin', 'ini', 'ini_roi', 'fin', 'fin_roi', ...
        'xx', 'yy');
end

clear F

%% run PIV and analyze results

% run piv
[xx, yy, uu, vv] = piv(ini, fin, ini_roi, fin_roi, xx, yy, samplen, ...
    sampspc, intrlen, npass, valid_max, valid_eps, lowess_span_pts, ...
    spline_tension, min_frac_data, min_frac_overlap, 1);

% compute strain values
[xgrid, ygrid] = meshgrid(xx, yy);
[displ, ~, ~, Dd, ~, ~, ~, ~] = ...
    deformation(xgrid, ygrid, uu, vv, ~isnan(uu));

% print and plot standard results
test_piv_util_plot(xx, yy, uu, vv, displ, Dd);
