% Run sandboxpiv on a series of images.  Parameters (detailed below) are
% expected to be saved in "input.mat" in the working folder. This funtion
% can be run from the Matlab path on any image series with locally defined
% parameters.
%
% SandboxPIV parameters ----------
%
% S = Vector. Width in pixels of sample window.  Each element corresponds
% to a single pass, # elements = # passes.
%
% SPC = Vector.  Distance between samples, # elements = # passes.
%
% ML = Vector.  Maximum displacement to the left (in pixels, in the image), # elements = # passes.
%
% MR = Vector.  " " right " "
%
% MU = Vector.  " " up " "
%
% MD = Vector.  " " down " "
%
% validate = Logical.  Turn vector valiation/interpolation on(1) or off(0)
%
% epsilon0 = Scalar.  Parameter to normalized median filter, can be
% omitted if validate=0;
%
% epsilonthresh = Scalar.  " "
%
% maskfrac = Scalar.  Fraction of sample window that must contain sand to
% proceed, otherwise mask
%
% Series parameters -------
% 
% imfilenames = Cell vector.  Each element contains the name of a file.
%
% intrapairstep = Integer.  Step between images in the same pair (1=adjacent) 
%
% interpairstep = Integer.  Step between image pairs (1=adjacent) 
%
% view = String. Select 'side' or 'top' view.
%
% savename = String.  Stub of filename to save results under.  The full name will be savename###.mat

% load parameters
load input.mat

params.S = S;
params.SPC = SPC;
params.ML = ML;
params.MR = MR;
params.MU = MU;
params.MD = MD;
params.CBC = CBC;
params.validate_nmed = validate_nmed;
params.epsilon0 = epsilon0;
params.epsilonthresh = epsilonthresh;
params.validate_nstd = validate_nstd;
params.nstd = nstd;
params.maskfrac = maskfrac;
params.view = view;
  
% loop: call PIV, save results
nim = numel(imfilenames);
count = 0;
nsteps = numel(1:interpairstep:nim-intrapairstep);
for i = 1:interpairstep:nim-intrapairstep

    count = count+1;
    params.IM = imfilenames([i,i+intrapairstep]);
    
    fprintf('\nRUNPIVSERIES: step %i of %i\n',count,nsteps);
    [x,y,u,v,cval,mask] = sandboxpiv(params);
    
    save(sprintf('%s%03i.mat',savename,count),'params','x','y','u','v','cval','mask');
    
end

% end
