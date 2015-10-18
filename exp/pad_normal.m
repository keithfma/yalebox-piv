% Script. Apply pad to match normal velocities outside the wedge to the wedge
% surface
%
% Pad is a function of the Euclidian distance to the boundary. If the boundary
% moves, this patterned pad will move also.

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

%% create pad: simple sinusoidal function of distance from the boundary
on = 1; 
if on 

    % parameters
    
    ampl = 0.4;
    period = 25; % pixels 
    
    % run
    pad_func = @(x) ampl*(sin(2*pi*x/period)+1)/2;
    
    ini_dist = bwdist(ini_mask);    
    ini_pad = ini;
    ini_pad(~ini_mask) = pad_func(ini_dist(~ini_mask));
    
    fin_dist = bwdist(fin_mask);    
    fin_pad = fin;
    fin_pad(~fin_mask) = pad_func(fin_dist(~fin_mask));
    
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

