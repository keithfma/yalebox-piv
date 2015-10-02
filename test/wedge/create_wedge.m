% Create wedge test input data using the image pair and parameters below

%% Define parameters 

% file names
ini_file = 'fault_ss_01_sidef_250.png';
fin_file = 'fault_ss_01_sidef_251.png';
out_file = 'fault_ss_01_sidef_250_251.mat';

% mask
mask_on = true;
hue_lim = [0, 0.5]; 
val_lim = [0, 0.5];
entr_lim = [0.5, 1];
entr_win = 9;
morph_rad = 15;

% equalization parameters
eql_on = false;
num_tiles = [100, 100]; 
clip_limit = 0.1;
num_bins = 1e4;

%% Create dataset

% load images
ini_rgb = imread(ini_file);  
fin_rgb = imread(fin_file);  
ini_hsv = rgb2hsv(ini_rgb);
fin_hsv = rgb2hsv(fin_rgb);
ini = ini_hsv(:,:,3);
fin = fin_hsv(:,:,3);

% create and apply masks, if enabled
if mask_on    
    mask_manual = yalebox_prep_mask_manual(ini_rgb);    

    ini_mask_auto = yalebox_prep_mask_auto(ini_hsv, ...
        hue_lim, val_lim, entr_lim, entr_win, morph_rad, true);
    fin_mask_auto = yalebox_prep_mask_auto(fin_hsv, ...
        hue_lim, val_lim, entr_lim, entr_win, morph_rad, true);    
    
    ini(~ini_mask_auto | ~mask_manual) = 0;
    fin(~fin_mask_auto | ~mask_manual) = 0;
end

% equalize histogram, if enabled
if eql_on   
    ini = adapthisteq(ini, 'NumTiles', num_tiles, 'NBins', num_bins, ...
            'ClipLimit', clip_limit, 'Range', 'full');
    fin = adapthisteq(fin, 'NumTiles', num_tiles, 'NBins', num_bins, ...
            'ClipLimit', clip_limit, 'Range', 'full');
end

% create fake coordinate vectors, in pixel coordinates
xx = 1:size(ini, 2);
yy = 1:size(ini, 1);

% save dataset
save(out_file, 'ini', 'fin', 'xx', 'yy');