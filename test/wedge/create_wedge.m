% Create wedge test input data using the image pair and parameters below

%% Define parameters 

% file names
dir_name = '~/Documents/dissertation/yalebox-exp-fault/data/fault_ss_01/image/crop_sidef/';
ini_file = 'fault_ss_01_sidef_250.png';
fin_file = 'fault_ss_01_sidef_251.png';
out_file = 'fault_ss_01_sidef_250_251.mat';

% mask
hue_lim = [0, 0.5]; 
val_lim = [0, 0.5];
entr_lim = [0.5, 1];
entr_win = 9;
morph_open_rad = 15;
morph_erode_rad = 10;

% intensity parameters
eql_nwin = 31;

%% Create dataset

% load images
ini_rgb = imread([dir_name, ini_file]);  
fin_rgb = imread([dir_name, fin_file]);  
ini_hsv = rgb2hsv(ini_rgb);
fin_hsv = rgb2hsv(fin_rgb);

% create masks   
mask_manual = yalebox_prep_mask_manual(ini_rgb);

ini_mask_auto = yalebox_prep_mask_auto(ini_hsv, ...
    hue_lim, val_lim, entr_lim, entr_win, morph_open_rad, morph_erode_rad, true);
fin_mask_auto = yalebox_prep_mask_auto(fin_hsv, ...
    hue_lim, val_lim, entr_lim, entr_win, morph_open_rad, morph_erode_rad, true);
    
ini_mask = ini_mask_auto & mask_manual;
fin_mask = fin_mask_auto & mask_manual;

% convert to intensity
ini = yalebox_prep_intensity(ini_hsv, ini_mask, eql_nwin, true);
fin = yalebox_prep_intensity(fin_hsv, fin_mask, eql_nwin, true);

% create fake coordinate vectors, in pixel coordinates
xx = 1:size(ini, 2);
yy = 1:size(ini, 1);

% save dataset
save(out_file, 'ini', 'fin', 'xx', 'yy');