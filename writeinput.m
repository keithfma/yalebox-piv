% Script to write the input variable required by runpivseries.m

% user defined parameters
S = [30];
SPC = [30]; 
ML = [40] ; 
MR = [40]; 
MU = [40]; 
MD = [40]; 
CBC.nnbr = [0]; % number of neighbors for CBC routine for each pass.  Valid options are 8,4,1,0
CBC.offsetfrac = 1/4; % offset in fraction of a sample window between center and neighbor windows.  Constant.
validate_nmed = 0; 
epsilon0 = 0.2;
epsilonthresh = 1;
validate_nstd = 1;
nstd = 6;
maskfrac = 0.33; 
intrapairstep = 1;
interpairstep = 1;
view = 'side';
savename = 'test_';
imfilenames = 'pRetroshear02SideF_*.png'; % supply stub here, the list is computed below 

% parse image filenames
imdir = dir(imfilenames);
imfilenames = cell(numel(imdir),1);
for i = 1:numel(imdir)
    imfilenames{i} = imdir(i).name;
end

% write input.mat
save input.mat S SPC ML MR MU MD CBC validate_nmed epsilon0 epsilonthresh validate_nstd nstd maskfrac intrapairstep interpairstep view savename imfilenames