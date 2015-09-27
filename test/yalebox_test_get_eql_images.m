% Apply adaptive histogram equalization to all the images

% define parameters
files = {'fault_ss_01_sidef_030.png', ...
        'fault_ss_01_sidef_031.png', ...
        'fault_ss_01_sidef_250.png', ...
        'fault_ss_01_sidef_251.png'};
num_tiles = [100, 100]; 
clip_limit = 0.1;
num_bins = 1e4;

% local histogram equalization (CLAHE) ignoring non-sand pixels

% create output file names
for i = 1:length(files)
    outfiles{i} = [files{i}(end-6:end-4), '_eql.mat'];
end
    
% setup environment
addpath('../');

for i = 1:length(files)
    
    % read image file
    rgb = imread(files{i});
    hsv = rgb2hsv(rgb);
    v = hsv(:,:,3);
    
    % equalize histogram
    eql = adapthisteq(v, ...
        'NumTiles', num_tiles, ...
        'NBins', num_bins, ...
        'ClipLimit', clip_limit, ...
        'Range', 'full');
     
    % clean up
    imshow(eql);
    pause
    close all
    
    % save results
    save(outfiles{i}, 'eql', 'num_tiles', 'num_bins', 'clip_limit');
    
end