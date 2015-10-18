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
