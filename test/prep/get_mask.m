% Create mask matrices (manual and auto) for the two test datasets.

% define parameters
files = {'fault_ss_01_sidef_030.png', ...
        'fault_ss_01_sidef_031.png', ...
        'fault_ss_01_sidef_250.png', ...
        'fault_ss_01_sidef_251.png'};
hue_lim = [0, 0.5]; 
val_lim = [0, 0.5];
entr_lim = [0.5, 1];
entr_win = 9;
morph_rad = 15;

% create output file names
for i = 1:length(files)
    outfiles{i} = [files{i}(end-6:end-4), '_mask.mat'];
end
    
% setup environment
addpath('../');

for i = 1:length(files)
    
    % read image file
    rgb = imread(files{i});
    hsv = rgb2hsv(rgb);

    % create manual mask, same for all images
    if i == 1
        mask_manual = yalebox_prep_mask_manual(rgb);
    end
    
    % generate automatic masks, different for each image
    mask_auto = yalebox_prep_mask_auto(hsv, hue_lim, val_lim, entr_lim, ...
        entr_win, morph_rad, true);
    
    % clean up
    pause
    close all
    
    % save results
    save(outfiles{i}, 'mask_manual', 'mask_auto', 'hue_lim', 'val_lim', ...
        'entr_lim', 'entr_win', 'morph_rad');   
    
end