function [] = datapostprocess(rawfilestub,maskfilestub,coordtransfile)
% function [] = datapostprocess(rawfilestub,maskfilestub,coordtransfile)
%
% Post-processing for raw sandboxpiv data:
%
% 1) Converts positions and displacement vectors to world coordinates 
%
% 2) Stores all step-specific data in individual MAT files (same name as
% inputs in a new folder (pwd/output/)
%
% 3) Stores the coordinate grids and transformation functions in a single
% MAT file (pwd/output/COORDS.mat)
%
% Arguments:
%
% rawfilestub = String.  Input to dir() that will yield a list of all
% files to be updated, for example 'RUN*.mat'.  Be careful that this will
% not include any files NOT to be processed.
%
% maskfilestub = String. INput to dir() that will yield a list of all full
% resolution mask files, as produced by YaleboxMask.m
%
% coordtransfile = String.  Filename of the output from GetWoco.
%
%
% Keith Ma, Yale University 2012

% preparations
slash = filesep; % system specific slash

if exist('output','dir') ~= 7 % create output directory if needed
    mkdir('output');
end

% load files
load(coordtransfile,'Jp2w','Jw2p','xw0','yw0');
rawfiles = dir(rawfilestub);
maskfiles = dir(maskfilestub);

if numel(rawfiles)~=(numel(maskfiles)-1) % break for bad inputs
    fprintf(2,'\nError: the number of data files returned by dir(rawfilestub) should be \none less than the number of mask files returened by dir(maskfilestub), \nbecause sandboxpiv operates on image pairs.\n\n');
    return
else
    nf = numel(rawfiles);
end

% save coordinates --------------------------------------------------

load(rawfiles(1).name,'x','y','params'); 

% convert xy to world
n = numel(x);
[nr nc] = size(x);
xvec = reshape(x,n,1);
yvec = reshape(y,n,1);
x = reshape(Jp2w(1,1).*xvec+Jp2w(1,2).*yvec+xw0, nr, nc);
y = reshape(Jp2w(2,1).*xvec+Jp2w(2,2).*yvec+yw0, nr, nc);

% flip the y-direction for later convenience
if strcmp(params.view,'side');
    y = flipud(y);
end

load(maskfiles(1).name,'mask'); 

% create world coordinate grids
n = numel(mask); 
[nr,nc] = size(mask);
[xPmask,yPmask] = meshgrid(1:nc,1:nr);
xvec = reshape(xPmask,n,1);
yvec = reshape(yPmask,n,1);
xMask = reshape(Jp2w(1,1).*xvec+Jp2w(1,2).*yvec+xw0, nr, nc);
yMask = reshape(Jp2w(2,1).*xvec+Jp2w(2,2).*yvec+yw0, nr, nc);

% create readme variable
COORDREADME.x = 'Double 2D matrix.  X-direction world coordinates for the displacement vectors in meters with the s-point at (0,0).  These positions are the center of the sample window used to compute the vector in sandboxpiv.m';
COORDREADME.y = 'Double 2D matrix.  Y-direction world coordinates for the displacement vectors in meters with the s-point at (0,0).  These positions are the center of the sample window used to compute the vector in sandboxpiv.m';
COORDREADME.xMask = 'Double 2D matrix.  X-direction world coordinates in meters for the original images with the s-point at (0,0).';
COORDREADME.yMask = 'Double 2D matrix.  Y-direction world coordinates in meters for the original images with the s-point at (0,0).';
COORDREADME.Jp2w = '2x2 matrix.  Jacobian coordinate transformation matrix for the conversion from pixel to world coordinates. [xW,yW] = Jp2w*[xP;yP]+[xw0;yw0]; (for coordinates).  [uW,vW] = Jp2w*[uP;vP]; (for displacements).';
COORDREADME.Jw2p = '2x2 matrix.  Jacobian coordinate transformation matrix for the conversion from world to pixel coordinates. [xP,yP] = Jw2p*([xP;yP]-[xw0;yw0]); (for coordinates).  [uP,vP] = Jw2p*[uP;vP]; (for displacements)';
COORDREADME.xw0 = 'Scalar.  Offset (translation) for conversion between pixel and world coordinates, see Jp2w and Jw2p for the form of the conversion.';
COORDREADME.yw0 = 'Scalar.  Offset (translation) for conversion between pixel and world coordinates, see Jp2w and Jw2p for the form of the conversion.';

save([pwd, slash, 'output', slash, 'COORDS.mat'],'x','y','xMask','yMask','Jp2w','Jw2p','xw0','yw0','COORDREADME');

% process each file -------------------------------------------------------


% create NEW readme variable
clear README
README.params = 'Struct. The parameters used in the sandboxpiv.m';
README.uX = 'Double 2D matrix.  X-direction component of the displacement vector in world coordinates (meters/step).';
README.uY = 'Double 2D matrix.  Y-direction component of the displacement vector in world coordinates (meters/step).';
README.cval = 'Double 2D matrix.  Value of the peak in the normalized cross correlation function used to identify displacements.';
README.iUXY = 'Logical 2D matrix. Mask array flagging sample windows for which vectors were computed (1) and sample windows which did not contain enough sand to compute a vector (0). Same size as xW,yW,uW,vW.' ;
README.iMask = 'Logical 2D matrix. Mask array flagging locations in the original image that contained sand (1) and locations the did not (0).  Same size as the original image.'; 
    
for i = 1:nf
    
    % process data ------------------------------------    
    fprintf('%s\t...\t',rawfiles(i).name);
    load(rawfiles(i).name,'x','y','u','v','cval','mask','params'); 
    
    n = numel(u);
    [nr nc] = size(u);
    uvec = reshape(u,n,1);
    vvec = reshape(v,n,1);
    uX = reshape(Jp2w(1,1).*uvec+Jp2w(1,2).*vvec, nr, nc);
    uY = reshape(Jp2w(2,1).*uvec+Jp2w(2,2).*vvec, nr, nc);
   
    % rename to new conventions
    iUXY = mask; 
    
    % flip y-direction for convenience
    if strcmp(params.view,'side');
        uX = flipud(uX);
        uY = flipud(uY);
        iUXY = flipud(iUXY);
        cval = flipud(cval);
    end
    
    % process masks --------------------------------------------------
    fprintf('%s\n',maskfiles(i).name);
    load(maskfiles(i).name); % loads mask    

    % rename to new conventions
    iMask = mask; 
    
    % save
    save([pwd, slash, 'output', slash, rawfiles(i).name],'params','uX','uY','cval','iUXY','iMask','README');
    
end



    