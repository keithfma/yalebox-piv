function [] = yalebox_piv_step()
% Skeleton for re-implementation of yalebox PIV analysis routine

% Arguments, input:
%
%   ini = 2D matrix, double, range 0 to 1, normalize grayscale image from
%       the start of the step to be analyzed.
%
%   fin = 2D matrix, double, range 0 to 1, normalize grayscale image from
%       the end of the step to be analyzed.
%
%   xx = Vector, double, increasing, x-direction coordinate vector, length
%       must match the columns in ini and fin.
%
%   yy = Vector, double, increasing, y-direction coordinate vector, length
%       must match the rows in ini and fin.
%
%   npass = Scalar, integer, number of PIV grid refinement passes
%
%   samplen = Vector, length === npass, double, side length of
%       the square sample window, elements must be odd so that windows can be
%       symmetric around the center point
%
%   xrez = Vector, length == npass, integer, grid points in the x-direction
%       for the output grid. If any element is set to 0, the number of points
%       will be chosen such that the aspect ratio is approximately 1 (in pixel
%       coordinates).
%
%   yrez = Vector, length == npass, integer, grid points in the x-direction
%       for the output grid. If any element is set to 0, the number of points
%       will be chosen such that the aspect ratio is approximately 1 (in pixel
%       coordinates).
%
%   umax = Vector, length == npass, maximum x-direction displacement
%       in world coordinates, used to set the size of the PIV search window.
%
%   umin = Vector, length == npass, minimum x-direction displacement
%       in world coordinates, used to set the size of the PIV search window, note that
%       this value will typically be negative to allow for displacements in the
%       negative x-direction.
%
%   vmax = Vector, length == npass maximum y-direction displacement
%       in world coordinates, used to set the size of the PIV search window.
%
%   vmin = Vector, length == npass minimum y-direction displacement
%       in world coordinates, used to set the size of the PIV search window, note that
%       this value will typically be negative to allow for displacements in the
%       negative y-direction.

%   
%   verbose = Scalar, logical, flag to enable (true) or disable (false) verbose
%       output messages.
%
% (NOTE: add a 'dryrun' mode that does everything but PIV, printing verbose outputs and returning ancillary variables)

%
% Arguments, output:
% 
% xx,yy,u,v
%
% References:

% validate = Scalar, logical, flag to enable (true) or disable (false) vector 
%   validation/interpolation
%
% eps_0 = Scalar, double, parameter to normalized median filter used for vector 
%   validation, ignored if validate == false;
%
% eps_thresh = Scalar, double  " "
%
% min_sand_frac = Scalar, double. Minimum fraction of sample window that must contain 
%   sand to proceed with PIV
%
% intra_pair_step = Integer.  Step between images in the same pair (1=adjacent) 
%
% inter_pair_step = Integer.  Step between image pairs (1=adjacent) 
%
% view = String. Select 'side' or 'top' view.

% debug { 
data_dir = '/home/kfm/Documents/dissertation/yalebox-exp-fault/data/fault_ss_01/piv/';
ini = flipud(double(imread([data_dir, 'img1.png']))/double(uint16(inf)));
fin = flipud(double(imread([data_dir, 'img2.png']))/double(uint16(inf)));
xx = fliplr((1:size(ini,2))/1e3);
yy = (1:size(ini,1))/1e3;
npass = 2;
samplen = [51, 25];
yrez = [100, 200];
xrez = [0, 0];
verbose = true;
umax = [20, 10];
umin = -umax+1;
vmax = [20, 10];
vmin = -vmax+1;
% } debug

print_input(verbose, 'input', ini, fin, xx, yy, npass, samplen, ...
    xrez, yrez, umin, umax, vmin, vmax); 

check_input(ini, fin, xx, yy, npass, samplen, xrez, yrez, umin, umax, vmin, ...
    vmax);

[xrez, yrez] = equalize_unknown_grid_dims(xrez, yrez, size(ini,1), size(ini,2));

[umin, umax] = uv_lim_world_to_pixel(umin, umax, xx);
[vmin, vmax] = uv_lim_world_to_pixel(vmin, vmax, yy);

print_input(verbose, 'preprocessed input', ini, fin, xx, yy, npass, ...
    samplen, xrez, yrez, umin, umax, vmin, vmax); 

end

function [] = check_input(ini, fin, xx, yy, npass, samplen, xrez, yrez, ...
                  umin, umax, vmin, vmax)
              
validateattributes(ini,...
    {'double'}, {'2d', 'real', 'nonnan', '>=', 0, '<=' 1}, ...
    mfilename, 'ini');
[nr, nc] = size(ini);
validateattributes(fin,...
    {'double'}, {'2d', 'real', 'nonnan', '>=', 0, '<=' 1, 'size', [nr, nc]}, ...
    mfilename, 'fin');
validateattributes(xx, ...
    {'double'}, {'vector', 'real', 'nonnan', 'numel', nc}, ...
    mfilename, 'xx');
validateattributes(yy, ...
    {'double'}, {'vector', 'real', 'nonnan', 'numel', nr}, ...
    mfilename, 'yy');
validateattributes(npass, ...
    {'numeric'}, {'scalar', 'integer', 'nonnegative'}, ...
    mfilename, 'npass');
validateattributes(samplen, ...
    {'numeric'}, {'numel', npass, 'integer', 'odd', 'positive', 'nonnan', }, ...
    mfilename, 'samplen');
validateattributes(xrez, ...
    {'numeric'}, {'numel', npass, 'integer', 'nonnegative', '<=', nc}, ...
    mfilename,  'xrez');
validateattributes(yrez, ...
    {'numeric'}, {'numel', npass, 'integer', 'nonnegative', '<=', nr}, ...
    mfilename, 'yrez');
validateattributes(umin, ...
    {'double'}, {'numel', npass}, ...
    mfilename, 'umin');
validateattributes(umax, ...
    {'double'}, {'numel', npass}, ...
    mfilename, 'umax');
validateattributes(vmin, ...
    {'double'}, {'numel', npass}, ...
    mfilename, 'vmin');
validateattributes(vmax, ...
    {'double'}, {'numel', npass}, ...
    mfilename, 'vmax');

end

function [xrez, yrez] = equalize_unknown_grid_dims(xrez, yrez, nr, nc)
% Approximately equalize any unknown grid dimensions 

for ii = 1:length(xrez)
    if xrez(ii) ~= 0 && yrez(ii) == 0
        % match unknown y-grid to known x-grid 
        pts = linspace(1, nc, xrez(ii));
        spc = pts(2)-pts(1);
        yrez(ii) = round(nr/spc)+1;        
    elseif xrez(ii) == 0 && yrez(ii) ~= 0
        % match unknown x-grid to known y-grid     
        pts = linspace(1, nr, yrez(ii));
        spc = pts(2)-pts(1);
        xrez(ii) = round(nc/spc)+1;    
    elseif xrez(ii) == 0 && yrez(ii) == 0
        % both grids are unknown, error    
        error('Unknown grid dimensions (both NaN) for pass %i', ii);
    end            
end

end

function [uvmin, uvmax] = uv_lim_world_to_pixel(uvmin, uvmax, xy)
% Convert displacement limits from world to pixel coordinates, one direction at
% a time. Sort to preserve the correct min/max regardless of the world
% coordinate axis polarity.

dxy = xy(2)-xy(1); 
uvminmax = sort([uvmin(:), uvmax(:)]/dxy, 2);
uvmin = uvminmax(:,1);
uvmax = uvminmax(:,2);

end

% verbose message subroutines --------------------------------------------------

function print_sep()
% Print a separator line for verbose output messages

fprintf('----------\n');

end

function print_input(verbose, msg, ini, fin, xx, yy, npass, samplen,...
                        xrez, yrez, umin, umax, vmin, vmax)
%
% Display values (or a summary of them) for the input arguments

if verbose
    print_sep;
    fprintf('%s\n', msg);
    fprintf('ini: size = [%i, %i], min = %.2f. max = %.2f, masked = %.2f%%\n',...
        size(ini, 1), size(ini, 2), min(ini(:)), max(ini(:)), ...
        sum(ini(:) ~= 0)/numel(ini)*100);
    fprintf('fin: size = [%i, %i], min = %.2f. max = %.2f, masked = %.2f%%\n',...
        size(fin, 1), size(fin, 2), min(fin(:)), max(fin(:)), ...
        sum(fin(:) ~= 0)/numel(fin)*100);
    fprintf('xx: length = %i, min = %.3f, max = %.3f, delta = %.3f\n', ...
        length(xx), min(xx), max(xx), xx(2)-xx(1));
    fprintf('yy: length = %i, min = %.3f, max = %.3f, delta = %.3f\n', ...
        length(yy), min(yy), max(yy), yy(2)-yy(1));
    fprintf('npass: %i\n', npass);
    fprintf('samplen: %s\n', sprintf('%i  ', samplen));
    fprintf('xrez: %s\n', sprintf('%i  ', xrez));
    fprintf('yrez: %s\n', sprintf('%i  ', yrez));
    fprintf('umin: %s\n', sprintf('%.2f  ', umin));
    fprintf('umax: %s\n', sprintf('%.2f  ', umax));
    fprintf('vmin: %s\n', sprintf('%.2f  ', vmin));
    fprintf('vmax: %s\n', sprintf('%.2f  ', vmax));
end

end

% % grid refinement loop (start)
% 
% % convert lengths from world to pixels
% 
% % compute sample window center points in pixels coords
% 
% % pad data and sample window center points (window+search_range+offsets, use padarray)
% 
% % (can I normalize away the penalty for flow out of the image?)
% 
% try 
%     %% Parse and check inputs, initialize
%     
%     % prompt user
%     fprintf('SANDBOXPIV VERSION %s\n',nversion);
%     
%     % interrogation variables
%     S = SETTINGS.S;
%     SPC = SETTINGS.SPC;
%     ML = SETTINGS.ML;
%     MR = SETTINGS.MR;
%     MU = SETTINGS.MU;
%     MD = SETTINGS.MD;
%     npass = numel(S); % number of grid refining passes
%     if min( [numel(SPC) numel(ML) numel(MR) numel(MU) numel(MD)]==npass )==0; error('Interrogation variable vectors must be the same length, exiting'); end
%     
%     % validation variables
%     validate_nmed = SETTINGS.validate_nmed;
%     validate_nstd = SETTINGS.validate_nstd;
%     
%     if validate_nmed==1
%         epsilon0 = SETTINGS.epsilon0;
%         epsilonthresh = SETTINGS.epsilonthresh;
%         
%         if numel(epsilon0)~=1 || numel(epsilonthresh)~=1; error('Bad value for .epsilon0 or .epsilonthresh, exiting'); end
%         
%     elseif validate_nmed==0
%         epsilon0 = NaN;
%         epsilonthresh = NaN;
%         
%     else error('Bad value for .validate_nmed, exiting\n');       
%     end
%     
%     if validate_nstd==1
%         nstd = SETTINGS.nstd;
%         
%         if numel(nstd)~=1; error('Bad value for .nstd, exiting'); end
%         
%     elseif validate_nstd==0
%         nstd = NaN;
%         
%     else error('Bad value for .validate_nstd, exiting');        
%     end
%     
%     % image variables
%     if numel(SETTINGS.IM)~=2 || ~ischar(SETTINGS.IM{1}) || ~ischar(SETTINGS.IM{1}); error('Bad value for .IM, exiting.'); end
%     
%     im1 = double(rgb2gray(imread(SETTINGS.IM{1})))/255; % assumes uint8!
%     im2 = double(rgb2gray(imread(SETTINGS.IM{2})))/255;
%     [nr nc] = size(im1);
%     
%     % Correlation-Based-Correction Variables
%     CBC = SETTINGS.CBC;    
%     for t = 1:numel(CBC.nnbr)
%         if max(CBC.nnbr(t)==[0 1 4 8])==0; error('Bad value for CBC.nnbr(%i). Exiting\n',t); end
%     end
%     
%     % other variables
%     maskfrac = SETTINGS.maskfrac;
%     if numel(maskfrac)~=1; error('Bad value for .maskfrac, exiting'); end
%     
%     if strcmp(SETTINGS.view,'side')==0 && strcmp(SETTINGS.view,'top')==0; error('Bad value for .view, exiting'); end
%     
%     
%     %% PIV
%     
%     % iterative grid refinement loop
%     for pass = 1:npass
%         
%         fprintf('-Pass %i of %i\n-window size = %ix%i, grid spacing = %i\n-correlation-based correction neighbors = %i, spacing = %i\n',...
%             pass,npass,S(pass),S(pass),SPC(pass),CBC.nnbr(pass),ceil(CBC.offsetfrac*S(pass)));
%         
%         CBCnwin = CBC.nnbr(pass)+1; % total number of windows to compute
%         
%         % define sampling grid
%         if pass>1 % mask previous output, store old positions for displacement interpolation
%             
%             u = u.*mask; % sets masked regions to 0
%             v = v.*mask;
%             oldsc = sc; % takes truncated windows into account
%             oldsr = sr;
%             oldu = u;
%             oldv = v;
%         end
%         
%         % TEST
%         % Align center of the first row of sample windows with the lower edge of the sand (upper edge of the image).
%         [sc sr] = meshgrid(1:SPC(pass):nc,nr:-SPC(pass):1); % sample window centers        
%         [nsr,nsc] = size(sr); % dimensions of sampling grid        
%         toEdgeFromCenter = round( (S(pass)-1) /2); % distance from the center of a sample window to the center rounded to a whole number 
%         sllc = sc(1,:)-toEdgeFromCenter; % lower left corners of sample cells in r,c coords, whole numbers
%         sllr = sr(:,1)-toEdgeFromCenter;         
%         
%         % initial guess for displacement field
%         if pass==1 % first pass, no information
%             u = zeros(nsr,nsc); % x-dir offset (intrinsic pixel units), ZEROS FOR NO DATA, SUBJECT TO CHANGE
%             v = zeros(nsr,nsc); % y-dir offset " "
%             
%         else % refine data from previous pass
%             
%             % crude interpolation to fill in outside the measured data (allows later use of linear,cubic,spline interpolation without introducing spurious edge gradients)
%             % must interpolate to the image boundaries as well
%             hasdata = oldu~=0 | oldv~=0; % use only points with measured data to construct the interpolant
%             olduTSI = TriScatteredInterp(oldsc(hasdata),oldsr(hasdata),oldu(hasdata),'nearest');
%             oldvTSI = TriScatteredInterp(oldsc(hasdata),oldsr(hasdata),oldv(hasdata),'nearest');
%             [intsc intsr] = meshgrid([1,sc(1,:),nc],[nr,sr(:,1)',1]); % old sample positions AND points on image boundaries, so no NaNs appear outside the convex hull of the data in later interpolation step
%             oldu = olduTSI(intsc,intsr);
%             oldv = oldvTSI(intsc,intsr);
%             
%             % improved interpolation to refined grid
%             u = interp2(intsc,intsr,oldu,sc,sr); % linear
%             v = interp2(intsc,intsr,oldv,sc,sr);
%             
%         end
%         
%         % setup correlation-based-correction (CBC)
%         CBCoffset = ceil(S(pass)*CBC.offsetfrac); % offset between central and adjacent sample windows for correlation plane correction
%         
%         switch CBC.nnbr(pass) % get offsets to each neighboring window
%             
%             case 8 % all 8 neighbors included
%                 [CBCoffr CBCoffc] = ind2sub([3,3],1:9);
%                 
%             case 4 % 4 neighbors included (NESW)
%                 [CBCoffr CBCoffc] = ind2sub([3,3],[2 4 5 6 8]);
%                 
%             case 1 % 1 neighbor included (E only, but this is arbitrary)
%                 [CBCoffr CBCoffc] = ind2sub([3,3],[5 8]);
%                 
%             case 0 % no CBR
%                 [CBCoffr CBCoffc] = ind2sub([3,3],5);
%                 
%         end
%         CBCoffr = (CBCoffr-2)*CBCoffset;
%         CBCoffc = (CBCoffc-2)*CBCoffset;
% 
%         % pad images with 0's to accomodate all possible displacements, take correlation plane summation into account
%         npadt = max( round( -min(v(:))+MD(pass)+CBCoffset+toEdgeFromCenter ), CBCoffset+toEdgeFromCenter ); % top pad = min vertical velocity + downwards search range + correlation plane offset
%         npadb = max( round( max(v(:))+MU(pass)+CBCoffset+toEdgeFromCenter ), CBCoffset+toEdgeFromCenter ); % bottom pad = max vertical velocity + upwards search range + " "
%         npadl = max( round( -min(u(:))+ML(pass)+CBCoffset+toEdgeFromCenter ), CBCoffset+toEdgeFromCenter ); % left pad = min leftward velocity + leftwards search range + " "
%         npadr = max( round( max(u(:))+MR(pass)+CBCoffset+toEdgeFromCenter ), CBCoffset+toEdgeFromCenter ); % right pad = max rightward velocity + rightwards search range + " "
%         
%         % ORIGINAL
% %         % pad images with 0's to accomodate all possible displacements, take correlation plane summation into account
% %         npadt = max( round( -min(v(:))+MD(pass)+CBCoffset ), CBCoffset ); % top pad = min vertical velocity + downwards search range + correlation plane offset
% %         npadb = max( round( max(v(:))+MU(pass)+CBCoffset ), CBCoffset ); % bottom pad = max vertical velocity + upwards search range + " "
% %         npadl = max( round( -min(u(:))+ML(pass)+CBCoffset ), CBCoffset ); % left pad = min leftward velocity + leftwards search range + " "
% %         npadr = max( round( max(u(:))+MR(pass)+CBCoffset ), CBCoffset ); % right pad = max rightward velocity + rightwards search range + " "
%         
%         im1pad = [zeros(npadt,nc+npadl+npadr);...
%             zeros(nr,npadl), im1, zeros(nr,npadr);...
%             zeros(npadb,nc+npadl+npadr)];
%         
%         im2pad = [zeros(npadt,nc+npadl+npadr);...
%             zeros(nr,npadl), im2, zeros(nr,npadr);...
%             zeros(npadb,nc+npadl+npadr)];
%                 
%         % allocate other vars
%         cval = zeros(nsr,nsc); % correlation value, not including subpixel offsets
%         mask = true(nsr,nsc); % one where displacement can be computed, 0 where it cannot
%         
%         % compute vectors for all samples
%         fprintf('--Computing displacements\n');
%         flag = false(nsr,nsc); % flags for discarded vectors
%         
%         
%         parfor i = 1:nsr
% %         for i = 1:nsr
%             for j = 1:nsc
%                 
%                 % fractional masking
%                 centerwinr = ( sllr(i) : sllr(i)+S(pass)-1 )+npadt; % sample window rows, [min max]
%                 centerwinc = ( sllc(j) : sllc(j)+S(pass)-1 )+npadl; % " " cols [min max]
%                 centerwin = im1pad(centerwinr,centerwinc);
%                 
%                 if sum(sum(centerwin~=0))/numel(centerwin) < maskfrac % don't bother computing if the central window doesn't contain sufficient sand
%                     mask(i,j) = 0;
%                     continue
%                 end
%                 
%                 % get offsets, assuming the full range exists - these are approximate for neighboring windows in the CBC routine (which is right)                
%                 offr = round(v(i,j))+(-MD(pass):MU(pass)); % Error corrected
%                 offc = round(u(i,j))+(-ML(pass):MR(pass));
%                 
%                 % get size of the correlation plane
%                 ncr = numel(offr);
%                 ncc = numel(offc);
%                 
%                 % reset vars
%                 localcor = zeros(ncr,ncc,9);
%                 
%                 % get correlation planes
%                 for p = 1:CBCnwin
%                     
%                     %  get sample window
%                     sampwinr = ( sllr(i) : sllr(i)+S(pass)-1 )+npadt+CBCoffr(p); % sample window rows, [min max]
%                     sampwinc = ( sllc(j) : sllc(j)+S(pass)-1 )+npadl+CBCoffc(p); % " " cols [min max]
%                     sampwin = im1pad(sampwinr,sampwinc);
%                     
%                     % fractional masking
%                     if sum(sum(sampwin~=0))/numel(sampwin) < maskfrac
%                         
%                         % continue % don't compute this window if it doesn't contain sufficient sand
%                                                
%                         % find a random position that DOES work, unless iteration becomes ridiculous
%                         randitr = 1;
%                         itrlimit = 1000;
%                         while randitr <= itrlimit 
%                             
%                             % choose a random position wihtin the allowed offset range
%                             randomoffsetr = round(rand(1)*CBCoffset);
%                             randomoffsetc = round(rand(1)*CBCoffset);
%                             
%                             %  get sample window
%                             sampwinr = ( sllr(i) : sllr(i)+S(pass)-1 )+npadt+randomoffsetr; % sample window rows, [min max]
%                             sampwinc = ( sllc(j) : sllc(j)+S(pass)-1 )+npadl+randomoffsetc; % " " cols [min max]
%                             sampwin = im1pad(sampwinr,sampwinc);
%                             
%                             if sum(sum(sampwin~=0))/numel(sampwin) >= maskfrac % success!
%                                 break
%                             else % keep looping 
%                                 randitr = randitr+1;
%                             end
%                             
%                             if randitr==itrlimit
%                                 fprintf('Skipping a CBC window, too many iteration trying to find a suitable window\n');
%                             end
%                             
%                             
%                         end
%                     
%                     end
%                     
%                     % get padded interrogation window
%                     intrwinr = (sampwinr(1)+min(offr):sampwinr(end)+max(offr));
%                     intrwinc = (sampwinc(1)+min(offc):sampwinc(end)+max(offc));
%                     intrwin = im2pad(intrwinr,intrwinc);
%                     
%                     fullcor = normxcorr2(sampwin,intrwin);
%                     localcor(:,:,p) = fullcor(S(pass):end-S(pass)+1,S(pass):end-S(pass)+1);
%                     
%                 end
%                                
%                 % add local correlation planes
%                 cor = sum(localcor,3);
%                 
%                 
%                 
%                 %%  Sub-pixel estimation
%                 
%                 % 9-point Gaussian (Nobach & Honkanen 2005, Experiments in Fluids)
%                 [~, ind] = max(cor(:)); % find max of corr plane
%                 [r,c] = ind2sub(size(cor),ind);
%                 
%                 if r~=1 && r~=ncr && c~=1 && c~=ncc % compute if the peak is in the interior of the correlation plane
%                     
%                     % offset the normalized correlation plane so there are no negative values as required by the gaussian peak model
%                     coroffset = abs(min(cor(:)));
%                     cor = cor+coroffset;
%                     
%                     % compute coefficients (could displose of loops, but it is easier this way!)
%                     c10 = 0; c01 = 0; c11 = 0; c20 = 0; c02 = 0; c00 = 0;
%                     for ii = -1:1
%                         for jj = -1:1
%                             logterm = log(cor(r+jj,c+ii));
%                             c10 = c10 + ii*logterm/6;
%                             c01 = c01 + jj*logterm/6;
%                             c11 = c11 + ii*jj*logterm/4;
%                             c20 = c20 + (3*ii^2-2)*logterm/6;
%                             c02 = c02 + (3*jj^2-2)*logterm/6;
%                             c00 = c00 + (5-3*ii^2-3*jj^2)*logterm/9;
%                         end
%                     end
%                     
%                     % compute sub-pixel displacement
%                     dr = ( c11*c10-2*c01*c20 )/( 4*c20*c02 - c11^2 );
%                     dc = ( c11*c01-2*c10*c02 )/( 4*c20*c02 - c11^2 );
%                     cval(i,j) = exp( c00-c20*dc^2-c11*dc*dr-c02*dr^2 )-coroffset; % remove offset from cval
%                     
%                     % apply subpixel displacement
%                     if abs(dr)<1 && abs(dc)<1 % subpixel estimation worked, there is a nice peak
%                         u(i,j) = offc(c)+dc;
%                         v(i,j) = offr(r)+dr;
%                     else % subpixel estimation failed, the peak is ugly and the displacement derived from it will stink
%                         flag(i,j) = 1;
%                         
%                     end
%                     
%                 else % drop vector and interpolate: lack of subpixel displacement will cause spurious gradients in the dataset
%                     flag(i,j) = 1;
%                     
%                 end
% 
% %                 % Whittaker-reconstruction (Lourenco & Krothapalli 1995,
% %                 % Experiments in Fluids) Thier implementation uses a 5*5
% %                 % region of the correlation plane to compute the
% %                 % interpolation at each point, and refines the grid
% %                 % iteratively.  Each iteration reduces the size of the grid
% %                 % by a factor of 2, so the current 6 iterations give a
% %                 % theoretical precision of 1/64th of the original window
% %                 % width.  It is a synch to use the whole correlation plane,
% %                 % which prevents the routine from failing, ever!
% %                 
% %                 % find max of corr plane
% %                 [~, ind] = max(cor(:));
% %                 [pkr,pkc] = ind2sub(size(cor),ind); % integer (pixel) position of the correlation peak
% %                 pkval = cor(pkr,pkc); % value of the correlation peak
% %                 
% %                 if pkr>=3 && pkr<=ncr-2 && pkc>=3 && pkc<=ncc-2 % compute if the peak is in the interior of the correlation plane
% %                     
% %                     % store initial integer position of the correlation peak
% %                     pixelpkr = pkr;
% %                     pixelpkc = pkc;
% %                     
% %                     % initialize whittaker interpolation routine
% %                     [corcols corrows] = meshgrid(1:ncc,1:ncr); % row,col coordinates of all points in the correlation plane
% %                     
% %                     % refine to at least 1e-3 pixels
% %                     nitr = ceil(log(1000*S(pass))/log(2));
% %                     
% %                     % use 5x5 grid
% %                     localvals = cor(pixelpkr-2:pixelpkr+2,pixelpkc-2:pixelpkc+2); localvals = localvals(:);
% %                     localcols = corcols(pixelpkr-2:pixelpkr+2,pixelpkc-2:pixelpkc+2); localcols = localcols(:);
% %                     localrows = corrows(pixelpkr-2:pixelpkr+2,pixelpkc-2:pixelpkc+2); localrows = localrows(:);
% %                     
% %                     % interpolate peak position, iteratively refining grid
% %                     for spitr = 1:nitr
% %                         
% %                         % find positions to interpolate to
% %                         interpgrdspc = 2^-spitr; % grid refines each iteration
% %                         [interpcorcols interpcorrows] = meshgrid([pkc-interpgrdspc,pkc,pkc+interpgrdspc],[pkr-interpgrdspc,pkr,pkr+interpgrdspc]); % row,col positions of the refined gridpoints to interpolate to
% %                         interpcor = nan(3); % declare
% %                         interpcor(2,2) = pkval; % the correlation value at the central point is known
% %                         
% %                         
% %                         % interpolate, whittaker style.
% %                         for k = [1:4,6:9] % all points excpet known central point
% %                             
% % %                             % use all values
% % %                             alpha = interpcorcols(k)-corcols(:);
% % %                             beta = interpcorrows(k)-corrows(:);
% % %                             interpcor(k) = sum( cor(:).*sinc(alpha).*sinc(beta) );
% %                             
% %                             % use 5x5 grid of values
% %                             alpha = interpcorcols(k)-localcols;
% %                             beta = interpcorrows(k)-localrows(:);
% %                             interpcor(k) = sum( localvals.*sinc(alpha).*sinc(beta) );
% %                             
% %                             
% %                         end
% %                         
% %                         % find refined peak position
% %                         [~, ind] = max(interpcor(:));
% %                         [intpkr,intpkc] = ind2sub(size(interpcor),ind); % integer (pixel) position of the correlation peak
% %                         pkval = interpcor(intpkr,intpkc); % value of the correlation peak
% %                         pkr = pkr+interpgrdspc*(intpkr-2);
% %                         pkc = pkc+interpgrdspc*(intpkc-2);
% %                         
% %                     end
% %                     
% %                     % apply subpixel displacement
% %                     u(i,j) = offc(pixelpkc)+(pkc-pixelpkc);
% %                     v(i,j) = offr(pixelpkr)+(pkr-pixelpkr);
% %                     
% %                 else % drop vector and interpolate: lack of subpixel displacement will cause spurious gradients in the dataset
% %                     flag(i,j) = 1;
% %                 end
%                 
%             end
%         end
%         
%         %% Vector Validation
%         
%         if validate_nmed || validate_nstd
%             
%             fprintf('--Vector validation\n');
%             
%             
%             if validate_nmed % SHOULD BE CHECK FOR BUGS
%                 
%                 fprintf('---Normalized median filter\n')
%                 
%                 % identify spurious vectors using normalized median filter of Westerweel & Scarano.
%                 for i = 1:nsr
%                     for j = 1:nsc
%                         
%                         if mask(i,j)==0 % skip if no vector
%                             continue
%                         end
%                         
%                         % extract values
%                         subu = u(max(i-1,1):min(i+1,nsr),max(j-1,1):min(j+1,nsc)); % center and 8-neighbors (if available)
%                         subv = v(max(i-1,1):min(i+1,nsr),max(j-1,1):min(j+1,nsc));
%                         center = subu==u(i,j);
%                         adju = subu(~center); adjv = subv(~center); % neighbors only
%                         hasdata = adju~=0 | adjv~=0; % drop any 0s
%                         adju = adju(hasdata); adjv = adjv(hasdata);
%                         
%                         if isempty(adju) || isempty(adjv); % no valid neighbors, skip
%                             continue
%                         end
%                         
%                         % get the median magnitude adjacent vector
%                         medind = ceil(numel(adju)/2); % choose the index of the median value (discrete) in a sorted list
%                         adjmag = sqrt( adju.^2+adjv.^2 );
%                         [~, sortind] = sort(adjmag);
%                         adju = adju(sortind); adjv = adjv(sortind);
%                         umed = adju(medind);
%                         vmed = adjv(medind);
%                         
%                         % get the median magnitude residual vector
%                         res = sqrt( (u(i,j)-adju).^2+(v(i,j)-adjv).^2 ); % magnitude of resisdual vectors
%                         [res] = sort(res);
%                         resmed = res(medind);
%                         
%                         testval = sqrt( (u(i,j)-umed)^2+(v(i,j)-vmed)^2 )/(resmed+epsilon0);
%                         
%                         if testval>epsilonthresh % bad vector
%                             flag(i,j) = 1;
%                         end
%                         
%                     end
%                 end
%                 
%             end
%             
%             if validate_nstd
%                 
%                 fprintf('---Standard deviation filter\n');
%                 
%                 meanu = mean(u(mask));
%                 meanv = mean(v(mask));
%                 stdu = std(u(mask));
%                 stdv = std(v(mask));
%                 flag(mask) = max(flag(mask), abs(u(mask)-meanu)>nstd*stdu | abs(v(mask)-meanv)>nstd*stdv);
%                 
%             end
%             
%         end
%         
%         % report percentage of vector discarded
%         fprintf('---Vectors discarded: %i, %.1f%%\n',sum(flag(:)),sum(flag(:))/sum(mask(:))*100);
%         
%         % drop rejected vectors
%         u(flag) = 0;
%         v(flag) = 0;
%         
%         
%     end
%     
%     %% Finishing steps
%     
%     
%     % remove masking at the base (s-point support, etc, these points should be filled via interpolation)
%     if strcmp(SETTINGS.view,'side')==1
%         mask = [ones(1,size(mask,2)); mask]; % paste a row of ones below (in world coords) the mask, blacked-out regions at the bed become holes
%         mask = bwfill(mask,'holes'); % fill the holes created above
%         mask = mask(2:end,:); % strip off the added row, holes are filled
%     end
%     
%     % interpolate dropped vectors ... step 1: extrapolate to pad data (build out its convex hull to contain all possible points)
%     pdist = SPC(pass); % distance out from data to create extrapolated pad
%     pr = [ (sr(1,:)+pdist)'; sr(:,end); (sr(end,:)-pdist)'; sr(:,1) ]; % positions of pad data points
%     pc = [ sc(1,:)'; sc(:,end)+pdist; sc(end,:)'; sc(:,1)-pdist ];
%     
%     hasdata = u~=0 | v~=0; % use only points with measured data to construct the interpolant
%     uTSI = TriScatteredInterp(sc(hasdata),sr(hasdata),u(hasdata),'nearest');
%     vTSI = TriScatteredInterp(sc(hasdata),sr(hasdata),v(hasdata),'nearest');
%     pu = uTSI(pc,pr); % extrapolate to pad
%     pv = vTSI(pc,pr);
%     
%     % interpolate dropped vectors ... step 2: interpolate using data and pad
%     uTSI = TriScatteredInterp([sc(hasdata); pc], [sr(hasdata); pr], [u(hasdata); pu],'linear');
%     vTSI = TriScatteredInterp([sc(hasdata); pc], [sr(hasdata); pr], [v(hasdata); pv],'linear');
%     u = uTSI(sc,sr);
%     v = vTSI(sc,sr);
%     
%     % mask displacement matrices
%     u = u.*mask;
%     v = v.*mask;
%     
%     % rename coordinate grids - the grids returned are the centers of the samping windows.
%     x = sc;
%     y = sr;
%    
%     
% catch err
%     disp(err)
%     keyboard
% end
% 
% % CUT
% 
% 
% %                     % replace masked areas with the mean (so these values go to 0 in the normalized correlation) -> SLIGHTLY WORSE RESULTS                  
% %                     sampwin(sampwin==0) = mean(sampwin(:));
% %                     intrwin(intrwin==0) = mean(intrwin(:));
%                     
% %                     % TOO SLOW - NEED A FASTER METHOD
% %                     % TEST: nearest neighbor fill for sample and interrogation windows
% %                     hasData = sampwin~=0;
% %                     [tempx,tempy] = meshgrid(1:size(sampwin,2),1:size(sampwin,1));
% %                     TSI = triScatteredInterp(tempx(hasData),tempy(hasData),sampwin(hasData),'nearest');
% %                     sampwin = TSI(tempx,tempy);
% %                     hasData = intrwin~=0;
% %                     [tempx,tempy] = meshgrid(1:size(intrwin,2),1:size(intrwin,1));
% %                     TSI = triScatteredInterp(tempx(hasData),tempy(hasData),intrwin(hasData),'nearest');
% %                     intrwin = TSI(tempx,tempy);
% 
% 
% 
% 
