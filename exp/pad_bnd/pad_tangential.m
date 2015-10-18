% Script. Apply pad to match tangental velocities outside the wedge with
% the surface velocities.
%
% Uses a nearest-neighbor interpolation to extend the intensity values
% along the surface of the wedge outwards along normal vectors. Since these
% edge pixels are translated along the wedge boundary, they constrain the
% tangental velocity.

%% load image pair data, get masks, pad edges

% parameters
input_file = 'fault_ss_01_sidef_250_251.mat';
npad = 50;

% run
load(input_file, 'ini', 'fin'); 

ini = padarray(ini, [npad, npad], 0, 'both');
fin = padarray(fin, [npad, npad], 0, 'both');

ini_mask = ini~=0;
fin_mask = fin~=0;

%% prepare images: adaptive histogram equalization first
% -> does not improve results, although the method is defective in that it does
% not do masked equalization

% parameters
tile_dim = 20;
num_tiles = floor(size(ini)/tile_dim);
nbins = 1000;

% run
ini = adapthisteq(ini, 'NumTiles', num_tiles, 'NBins', nbins, 'Range', 'full');
ini(~ini_mask) = 0;

fin = adapthisteq(fin, 'NumTiles', num_tiles, 'NBins', nbins, 'Range', 'full');
fin(~fin_mask) = 0;

%% create pad: simple nearest neighbor
% -> too much noise on the boundary, will not be useful
on = 0; 
if on 
    fprintf('simple nearest neighbor\n');

    % run
    [ini_dist, ini_nearest] = bwdist(ini_mask);
    ini_pad = ini(ini_nearest);
    
    [fin_dist, fin_nearest] = bwdist(fin_mask);
    fin_pad = fin(fin_nearest);
    
end
%% create pad: nearest neighbor with a boundary values taken from a surface "skin depth"
% -> values are more muted, still appears too noisy to be useful
on = 1;
if on
    fprintf('nearest neighbor with skin depth\n');
    
    % parameters    
    nskin = 20;
    
    % run
    ini_pad = ini;  
    fin_pad = fin;
    
    for ii = 0:nskin-1
        
        ini_mask_erode = imerode(ini_mask, strel('disk', ii));
        [~, ini_nearest_erode] = bwdist(ini_mask_erode);
        ini_pad = ini_pad+ini(ini_nearest_erode);        
       
        fin_mask_erode = imerode(fin_mask, strel('disk', ii));
        [~, fin_nearest_erode] = bwdist(fin_mask_erode);
        fin_pad = fin_pad+fin(fin_nearest_erode);               
        
    end
    ini_pad = ini_pad/nskin;
    fin_pad = fin_pad/nskin;


end
%% visualize results

% parameters
clim = [0, 0.4];

% run
figure;
while 1
    imagesc(ini_pad);
    caxis(clim);
    axis equal
    axis tight
    pause(1);    
    imagesc(fin_pad);
    caxis(clim);
    axis equal
    axis tight
    pause(1);
end