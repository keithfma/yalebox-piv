function [ini, fin, ini_mask, fin_mask, xx, yy] = ...
    piv_test_create_wedge(dir_name, ini_file, fin_file, hue_lim, val_lim, ...
        entr_lim, entr_win, morph_open_rad, morph_erode_rad, eql_nwin)
% function [ini, fin, ini_mask, fin_mask, xx, yy] = ...
%     piv_test_create_wedge(dir_name, ini_file, fin_file, hue_lim, val_lim, ...
%         entr_lim, entr_win, morph_open_rad, morph_erode_rad, eql_nwin)
%
% Generate input variables for PIV test run using experimental wedge images.
% This function performs the usual pre-processing steps with the given (or
% default) parameters.
%
% %

%% initialize

% set defaults
narginchk(3, 10);
if nargin < 4 || isempty(hue_lim)
   hue_lim = [0, 0.5]; 
end
if nargin < 5 || isempty(val_lim)
    val_lim = [0, 0.5];
end
if nargin < 6 || isempty(entr_lim)
    entr_lim = [0.5, 1];
end
if nargin < 7 || isempty(entr_win)
    entr_win = 9;
end
if nargin < 8 || isempty(morph_open_rad)
    morph_open_rad = 15;
end
if nargin < 9 || isempty(morph_erode_rad)
    morph_erode_rad = 10;
end
if nargin < 10 || isempty(eql_nwin)
    eql_nwin = 31;
end

% check for sane inputs (many checks are handled by subfunctions
assert(exist(dir_name, 'dir')==7, 'Directory %s does not exist\n', dir_name);
ini_file = [dir_name, filesep, ini_file];
fin_file = [dir_name, filesep, fin_file];
assert(exist(ini_file, 'file')==2, 'File %s does not exist\n', ini_file);
assert(exist(fin_file, 'file')==2, 'File %s does not exist\n', fin_file);

%% prepare images

% load images
ini_rgb = imread(ini_file);  
fin_rgb = imread(fin_file);  
ini_hsv = rgb2hsv(ini_rgb);
fin_hsv = rgb2hsv(fin_rgb);
ini = ini_hsv(:,:,3);
fin = fin_hsv(:,:,3);

% create masks   
mask_manual = prep_mask_manual(ini_rgb);

ini_mask_auto = yalebox_prep_mask_auto(ini_hsv, ...
    hue_lim, val_lim, entr_lim, entr_win, morph_open_rad, morph_erode_rad, 0);
fin_mask_auto = yalebox_prep_mask_auto(fin_hsv, ...
    hue_lim, val_lim, entr_lim, entr_win, morph_open_rad, morph_erode_rad, 0);
    
ini_mask = ini_mask_auto & mask_manual;
fin_mask = fin_mask_auto & mask_manual;

% convert to intensity
ini = prep_intensity(ini, ini_mask, eql_nwin, 0);
fin = prep_intensity(fin, fin_mask, eql_nwin, 0);

% create fake coordinate vectors, in pixel coordinates
xx = 1:size(ini, 2);
yy = 1:size(ini, 1);
