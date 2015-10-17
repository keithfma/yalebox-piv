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

%% create pad: simple nearest neighbor
% -> too much noise on the boundary, will not be useful
on = 1; 
if on 

    % run
    [ini_dist, ini_nearest] = bwdist(ini_mask);
    ini_pad = ini(nearest(ini_nearest));
    
    [fin_dist, fin_nearest] = bwdist(fin_mask);
    fin_pad = fin(nearest(fin_nearest));
    
end
%% create pad: nearest neighbor with a skin depth
%
on = 0;
if on
    
    
    skin_depth = 5;
    
    ini_pad_0 = ini;
    
    for ii = 0:skin_depth-1
        
        
        
        
    end


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