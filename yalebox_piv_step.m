function [xx, yy, uu, vv] = yalebox_piv_step()
% Re-implementation of yalebox PIV analysis routine
%
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
%   samplen = Vector, length === npass, integer, side length of the square
%       sample window
%
%   sampspc = Vector, length === npass, integer, spacing between adjacent sample
%       points in the (square) sample grid for each pass
%
%   umax, umin = Vector, length == npass, maximum x-direction displacement
%       in world coordinates, used to set the size of the PIV search window.
%       Values may be negative to allow for displacements in the negative
%       x-direction.
%
%   vmax, vmin = Vector, length == npass maximum y-direction displacement
%       in world coordinates,  used to set the size of the PIV search window.
%       Values may be negative to allow for displacements in the negative
%       x-direction.
%
%   verbose = Scalar, integer, flag to set verbosity level, (0) no verbose
%       output, (1) enable verbose text output messages, (2) enable debugging
%       plots
%
% Arguments, output:
%
%   xx, yy = Vector, double, coordinate vectors for the final output sample
%       grid, in world coordinate units
%
%   uu, vv = 2D matrix, double, computed displacement in the x- and y-directions
%       in world coordinate units
% 
%   smoothing factors?
%   some measure of quality?
%
% References:
%
% [2] Nobach, H., & Honkanen, M. (2005). Two-dimensional Gaussian
% regression for sub-pixel displacement estimation in particle image
% velocimetry or particle position estimation in particle tracking
% velocimetry. Experiments in Fluids, 38(4), 511-515.
% doi:10.1007/s00348-005-0942-3
%
% [3] Garcia, D. (2010). A fast all-in-one method for
% automated post-processing of PIV data. Experiments in Fluids, 50(5),
% 1247?1259. doi:10.1007/s00348-010-0985-y
%
% [4] Garcia, D. (2010). Robust smoothing of gridded data in one and higher
% dimensions with missing values. Computational Statistics & Data Analysis,
% 56(6), 2182. doi:10.1016/j.csda.2011.12.001
 
% NOTE: use NaNs instead of zeros to indicate the mask/roi

% debug {

% % single pass, test 01
% load('test01_input.mat', 'ini', 'fin', 'xx', 'yy');
% npass = 1;
% samplen = 30;
% sampspc = 15;
% umax =  0.05;
% umin = -0.05;
% vmax =  0.05;
% vmin = -0.05;
% ncbc = 9;
% verbose = 2;

% % single pass, test 02
% load('test02_input.mat', 'ini', 'fin', 'xx', 'yy');
% npass = 1;
% samplen = 50;
% sampspc = 25;
% umax =  0.02;
% umin = -0.02;
% vmax =  0.02;
% vmin = -0.02;
% ncbc = 9;
% verbose = 1;

% dual pass, test 01
load('test01_input.mat', 'ini', 'fin', 'xx', 'yy');
npass = 2;
samplen = [30, 20];
sampspc = [15, 10];
umax = [ 0.05,  0.03];
umin = [-0.05, -0.03];
vmax = [ 0.05,  0.03];
vmin = [-0.05, -0.03];
ncbc = [9, 9];
verbose = 1;

% % dual pass, test 02
% load('test02_input.mat', 'ini', 'fin', 'xx', 'yy');
% npass = 2;
% samplen = [50, 30];
% sampspc = [25, 15];
% umax = [ 0.02,  0.005];
% umin = [-0.02, -0.005];
% vmax = [ 0.02,  0.005];
% vmin = [-0.02, -0.005];
% ncbc = [9, 9];
% verbose = 1;

% % tri pass, test 02
% load('test02_input.mat', 'ini', 'fin', 'xx', 'yy');
% npass = 3;
% samplen = [60, 40, 20];
% sampspc = [30, 20, 10];
% umax = [ 0.02,  0.005,  0.0025];
% umin = [-0.02, -0.005, -0.0025];
% vmax = [ 0.02,  0.005,  0.0025];
% vmin = [-0.02, -0.005, -0.0025];
% ncbc = [9, 9, 9];
% verbose = 1;

% } debug

print_sep('INPUT ARGUMENTS', verbose);
print_input(ini, fin, xx, yy, npass, samplen, sampspc, umin, umax, ...
        vmin, vmax, ncbc, verbose);

check_input(ini, fin, xx, yy, npass, samplen, sampspc, umin, umax, ...
    vmin, vmax, ncbc);

[umin, umax] = uv_input_world_to_pixel(umin, umax, xx);
[vmin, vmax] = uv_input_world_to_pixel(vmin, vmax, yy);

[nr0, nc0] = size(ini);
[rr, cc] = sample_grid(sampspc(1), nr0, nc0);
uu = zeros(length(rr), length(cc));
vv = zeros(length(rr), length(cc));

% loop over PIV passes
for pp = 1:npass
    
    print_sep(sprintf('PIV pass %i of %i', pp, npass), verbose);
    print_pass(rr, cc, umax(pp), umin(pp), vmax(pp), vmin(pp), ncbc(pp), ...
        verbose);
    
    % init per-pass variables
    nr = length(rr);
    nc = length(cc);
    nv = vmax(pp)-vmin(pp);
    nu = umax(pp)-umin(pp);    
    
    roi = true(nr, nc);
    xcr_stack = nan(nv, nu, nr*nc);
    vorigin_stack = nan(nr*nc);
    uorigin_stack = nan(nr*nc);
                                          
    % loop over sample grid, gather cross-correlation data   
    for jj = 1:nc
        for ii = 1:nr
            
            % get linear index for ii, jj
            kk = ii+(jj-1)*nr;
            
            % skip if window center lies outside the roi at start or finish
            rchk = round(rr(ii));
            cchk = round(cc(jj));            
            if ini(rchk, cchk) == 0 || fin(rchk,cchk) == 0
                roi(ii, jj) = false;
                continue
            end
            
            % get sample and interrogation windows
            [samp, spos, intr, ipos, vorigin_stack(kk), uorigin_stack(kk)] = ...
                get_win(ini, fin, rr(ii), cc(jj), samplen(pp), ...
                    vv(ii,jj), vmin(pp), vmax(pp), ...
                    uu(ii,jj), umin(pp), umax(pp));
                
            % skip if interrogation window is empty
            if max(intr(:)) == 0
                roi(ii, jj) = false;
                continue
            end
                
            % compute correlation, trimming to valid range (see help)            
            xcr_stack(:, :, kk) = get_cross_corr(samp, intr);
           
            % plot_sample_point(ini, fin, samp, samppos, intr, intrpos, ...
            %     xcr_stack(:, :, kk), rpeak, cpeak, verbose);
           
        end
    end
    
    % loop over sample grid, compute displacements
    for jj = 1:nc
        for ii = 1:nr
            
            % get linear index for ii, jj
            kk = ii+(jj-1)*nr;
            
            % skip if indicated
            if roi(ii, jj) == false
                continue
            end
            
            % get subscripts for local correlation planes to include in analysis            
            [ixcr, jxcr] = get_cbc_stencil(ii, jj, ncbc(pp));
            
            % delete non-existant points
            valid = ixcr >= 1 & ixcr <= nr & jxcr >= 1 & jxcr <= nc;
            ixcr = ixcr(valid);
            jxcr = jxcr(valid);
            
            % convert to linear indices
            kxcr = ixcr+(jxcr-1)*nr;
            
            % delete points outside the roi
            valid = roi(kxcr);
            kxcr = kxcr(valid);
                        
            % combine correlation planes
            % multiply
            xcr = ones(nv, nu);
            for i = 1:length(kxcr)
                xcr = xcr.*xcr_stack(:, :, kxcr(i));
            end
                        
            % find the correlation plane maximum with subpixel accuracy
            
            [rpeak, cpeak, status] = find_peak(xcr_stack(:, :, kk));
            if status == false
                vv(ii, jj) = NaN;
                uu(ii, jj) = NaN;
                continue
            end     
     
            % get displacement in pixel coordinates
            vv(ii, jj) = vorigin_stack(kk)+rpeak-1;
            uu(ii, jj) = uorigin_stack(kk)+cpeak-1;
            
        end
    end
    
    % post-process and prep for next pass
    pp_next = min(npass, pp+1); % last pass uses same grid
    [rr_next, cc_next] = sample_grid(sampspc(pp_next), nr0, nc0);            
    [uu, vv] = post_process(uu, vv, roi, rr, cc, rr_next, cc_next);
    rr = rr_next;
    cc = cc_next;
    
end

% set interpolated values outside the roi to NaN
uu(~roi) = NaN;
vv(~roi) = NaN;

% convert displacements to world coordinates
uu = uu.*(xx(2)-xx(1));
vv = vv.*(yy(2)-yy(1));

% get world coordinate vectors for final sample grid
yy = interp1(1:nr0, yy, rr);
xx = interp1(1:nc0, xx, cc);

end

function check_input(ini, fin, xx, yy, npass, samplen, sampspc, umin, ...
             umax, vmin, vmax, ncbc)
% Check for sane input argument properties, exit with error if they do not
% match expectations.
              
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
    {'numeric'}, {'numel', npass, 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'samplen');
validateattributes(sampspc, ...
    {'numeric'}, {'numel', npass, 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'sampspc');
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
validateattributes(ncbc, ...
    {'double'}, {'numel', npass}, ...
    mfilename, 'ncbc');
for i = 1:npass
    assert(ismember(ncbc(i), [1, 3.1, 3.2, 5, 9]), ...
        'Invalid value of %g for ncbc(%i), options are 1, 3.1, 3.2, 5, 9', ...
        ncbc(i), i);
end

end

function [uvmin, uvmax] = ...
	uv_input_world_to_pixel(uvmin, uvmax, xy)
% Convert displacement limits from world to pixel coordinates, one
% direction at a time. Sort to preserve the correct min/max regardless of
% the world coordinate axis polarity.
%
% Arguments:
%
%   uvmin, uvmax = Scalar, integer, minimum and maximum displacement in
%       world coordinates for either x- or y-direction, rounded to integer
%       values away from zero.
%
%   xy = Vector, double, world coordinate vector for either x- or
%       y-direction

dxy = xy(2)-xy(1); 

uvminmax = sort([uvmin(:), uvmax(:)]/dxy, 2);
uvmin = floor(uvminmax(:,1));
uvmax = ceil(uvminmax(:,2));

end

function [rr, cc] = sample_grid(spc, nr0, nc0)
% Create sample grid using equally spaced points at integer pixel
% coordinates, with the remainder split evenly between edges.
%
% Arguments:
%
%   spc = Scalar, integer, grid spacing in pixels
%
%   nro, nc0 = Scalar, integer, grid dimensions (rows, columns) for the
%       original input data matrices (e.g. ini and fin).
%
%   rr, cc = Vector, integer, coordinate vectors for the sample grid in the
%       y/row and x/column directions, in pixels

rrem = mod(nr0-1, spc);
rri = 1+floor(rrem/2);
rr = rri:spc:nr0;

crem = mod(nc0-1, spc);
cci = 1+floor(crem/2);
cc = cci:spc:nc0;

end

function [swin, spos, iwin, ipos, vorigin, uorigin] = ...
    get_win(sdata, idata, rpt, cpt, slen, v0, v1min, v1max, u0, u1min, u1max)
% Extract the sample and interrogation windows for a given point. Returns the
% windows (padded as needed) and the displacement in pixel coordinates of the
% origin (element 1,1) of the valid correlation matrix of swin and iwin. The
% latter are used to convert peak position in the correlation plane to
% displacements in pixel coordiantes.
%
% Arguments:
% 
% sdata, idata = 2D matrix, double, initial and final data matrices from which
%   to extract the sample and interrogation windows (respectively).
%
% rpt, cpt = Scalar, double, location of the sample window centerpoint in pixel
%   coordinates.
%
% slen = Scalar, integer, number of points along a side of the square sample
%   window.
%
% u0, v0 = Scalar, double, estimated displacements from prior pass or initial
%   guess, defines the location of the interrogation window
%
% u1min, u1max, v1min, v1max = Scalar, integer, maximum and minimum
%   displacements in addition to the estimated displacements (v0, u0) in both
%   directions in pixel coordinates, defines the size and location of the
%   interrogation window.
%
% swin, iwin = 2D matrix, double, the sample and interrogation windows
%
% spos, ipos = Vector, double, position vectors for the sample and interrogation
%   windows, of the form [left, bottom, width, height], used by verbose plotting
%   routine
% 
% uorigin, vorigin = displacement in pixel coordinates of the origin (element
%   1,1) of the correlation matrix of swin and iwin, assumes size is the 'valid'
%   extent (a la conv2)

% parameters
hwidth = (slen-1)/2;

% get sample window index range, shifted to nearest whole pixel
rmin = rpt-hwidth;
rmax = rpt+hwidth;
radj = round(rmin)-rmin;
rmin = rmin+radj;
rmax = rmax+radj;

cmin = cpt-hwidth;
cmax = cpt+hwidth;
cadj = round(cmin)-cmin;
cmin = cmin+cadj;
cmax = cmax+cadj;

spos = [cmin, rmin, cmax-cmin, rmax-rmin];

% extract sample window, including pad if needed
swin = get_padded_subset(sdata, rmin, rmax, cmin, cmax);

% get interrogation window index range, shifted down/left to whole pixel
rmin = rpt+v0+v1min-hwidth;
rmax = rpt+v0+v1max+hwidth;
radj = floor(rmin)-rmin;
rmin = round(rmin+radj); % round to deal with floating point imprecision
rmax = round(rmax+radj);

cmin = cpt+u0+u1min-hwidth;
cmax = cpt+u0+u1max+hwidth;
cadj = floor(cmin)-cmin;
cmin = round(cmin+cadj); % round to deal with floating point imprecision
cmax = round(cmax+cadj);

ipos = [cmin, rmin, cmax-cmin, rmax-rmin];

% extract interrogation window, including pad if needed
iwin = get_padded_subset(idata, rmin, rmax, cmin, cmax);

% get displacement origin
vorigin = rmin-rpt+floor(slen/2);
uorigin = cmin-cpt+floor(slen/2);

end

function [win] = get_padded_subset(data, r0, r1, c0, c1)
% Extract subset of data in the range (r0:r1, c0:c1), padding with zeros
% where the indices extend beyond the limits of the data
%
%   data = 2D Matrix, data from which a subset is to be extracted
%
%   r0, r1 = Scalar, integer, rows requested for the subset, if these lie
%       outside the range [1, size(data,1], the output matrix will be padded
%       with zeros to maintain the requested size
%
%   c0, c1 = Scalar, integer, columns requested for the subset, if these lie
%       outside the range [1, size(data,1], the output matrix will be padded
%       with zeros to maintain the requested size

% get pad size and restrict window indices to valid range
pl = max(0, 1-c0);
c0 = max(1, c0);

nc = size(data, 2);
pr = max(0, c1-nc);
c1 = min(nc, c1);

pb = max(0, 1-r0);
r0 = max(1, r0);

nr = size(data, 1);
pt = max(0, r1-nr);
r1 = min(nr, r1);

% extract data and add pad
sub = data(r0:r1, c0:c1);
[snr, snc] = size(sub);
win = [zeros(pb, pl+snc+pr);
       zeros(snr, pl), sub, zeros(snr, pr);
       zeros(pt, pl+snc+pr)];  
   
% % debug {
% imagesc(data); caxis([0,1]); axis equal; hold on;
% plot([c0, c0, c1, c1, c0], [r0, r1, r1, r0, r0], 'LineWidth', 2, 'Color', 'k');
% hold off; drawnow; pause(0.01);
% % } debug
    
end

function [xcorr] = get_cross_corr(aa, bb)
% Compute normalized cross correlation, and crop to the 'valid' extent of the
% correlation (a la conv2). Allows for non-square aa, although that is not
% needed at this time.
%
% Arguments:
%   aa = 2D matrix, double, smaller 'template' matrix, as used here, this
%       is the sample window
%
%   bb = 2D matrix, double, larger matrix, as used here, this is the
%       interrogation window

fullxcorr = normxcorr2(aa, bb);

% compute pad size in both dimensions
aaSize = size(aa);
npre = aaSize;
npost = aaSize;

xcorr = fullxcorr( (1+npre(1)):(end-npost(1)+1), (1+npre(2)):(end-npost(2))+1);
            
end

function [rind, cind] = get_cbc_stencil(ii, jj, name)
% get subscripts for the neighboring sample points, which may lie outside the
% sample grid
%
% Arguments:
%
%   ii, jj = Scalar, integer, row and column position of the center point 
%
%   name = Scalar, double, numeric flag seleting the CBC stencil (e.g. 1 for no
%       CBC, 3.2 for CBC with 3 points along the 2nd dimensions (columns)

switch name
    case 1
        % no CBC
        rind = ii;
        cind = jj;
        
    case 3.1
        % 3-point stencil, along rows
        rind = [ii, ii  , ii  ];
        cind = [jj, jj-1, jj+1];
    
    case 3.2
        % 3-point stencil, along columns
        rind = [ii, ii-1, ii+1];
        cind = [jj, jj  , jj  ];
        
    case 5        
        % 5-point stencil
        rind = [ii, ii  , ii  , ii-1, ii+1];
        cind = [jj, jj-1, jj+1, jj  , jj  ];
        
    case 9 
        % 9-point stencil
        rind = [ii-1, ii  , ii+1, ii-1, ii  , ii+1, ii-1, ii  , ii+1];
        cind = [jj+1, jj+1, jj+1, jj  , jj  , jj  , jj-1, jj-1, jj-1];
end

end

function [rpk, cpk, stat] = find_peak(zz)
% Find the position of the peak in matrix zz with subpixel accuracy. Peakl
% location is determined from an explicit solution of two-dimensional
% Gaussian regression (see [2]). If the peak cannot be fit at subpixel
% accuracy, no peak is returned (see Arguments). This choice reflects the
% fact that a lack of subpixel displacement causes spurious gradients - it
% is preferable to drop the vector and interpolate.
%
% Arguments:
%   zz = 2D matrix, data plane in which to locate the peak
%
%   rpk = Scalar, double, row-coordinate location of the peak, set to -1 if
%       the peak cannot be fit.
%
%   cpk = Scalar, double, column-coordinate location of the peak, set to -1
%       if the peak cannot be fit
%
%   stat = Logical, scalar, return status flag, true if successful, false if
%       unsuccessful

[rpk, cpk] = find(zz == max(zz(:)));

% check failure conditions
%   1-2) no unique maximum
%   3-6) peak at the edge of the matrix
stat = true;
if numel(rpk) ~= 1 || numel(cpk) ~= 1 ...
        || rpk == 1 || rpk == size(zz, 1) || cpk == 1 || cpk == size(zz,2)        
    stat = false;
    return
end
    
% offset to eliminate non-positive (gaussian is always positive)
zz = zz-min(zz(:))+eps;

% compute coefficients 
c10 = 0; 
c01 = 0; 
c11 = 0; 
c20 = 0; 
c02 = 0; 
c00 = 0;
for ii = -1:1
    for jj = -1:1
        logterm = log(zz(rpk+jj,cpk+ii));
        c10 = c10 + ii*logterm/6;
        c01 = c01 + jj*logterm/6;
        c11 = c11 + ii*jj*logterm/4;
        c20 = c20 + (3*ii^2-2)*logterm/6;
        c02 = c02 + (3*jj^2-2)*logterm/6;
        c00 = c00 + (5-3*ii^2-3*jj^2)*logterm/9;
    end
end
                     
% compute sub-pixel displacement
dr = ( c11*c10-2*c01*c20 )/( 4*c20*c02 - c11^2 );
dc = ( c11*c01-2*c10*c02 )/( 4*c20*c02 - c11^2 );

% apply subpixel displacement
if abs(dr) < 1 && abs(dc) < 1
    % subpixel estimation worked, there is a nice peak
    rpk = rpk+dr;
    cpk = cpk+dc;
    
else
    % subpixel estimation failed, the peak is ugly and the displacement derived from it will stink
    stat = false;
end

end

function [uu1, vv1, sf] = post_process(uu0, vv0, roi0, rr0, cc0, rr1, cc1)
% Post-process PIV data using DCT-PLS to validate, replace and smooth
% vectors, and interpolate the results to the new sample grid. The
% post-processing algorithm is also used to extrapolate the data so that
% the current pass can be mapped to the new grid. The no-data region
% outside the roi is set to NaN, and a 1-point pad of NaNs are added prior
% to post-processing in order to accomplish the extrapolation.
% 
% Arguments:
%
%   uu0, vv0 = 2D matrix, double, displacement components for the current
%       sample grid in pixel coordinates.
%
%   roi0 = 2D matrix, logical, mask matrix for the current sample grid
%       indicating where there is data (1) and where there is none (0)
%
%   rr0, cc0 = Vector, integer, coordinate vectors for the current sample
%       grid y/row and x/colum dimensions, in pixels
%
%   rr1, cc1 = Vector, integer, coordinate vectors for the new sample grid 
%       y/row and x/colum dimensions, in pixels
%
%   uu1, vv1 = 2D matrix, double, displacement components for the new
%       sample grid in pixel coordinates.
%
%   sf = Scalar, double, smoothing factor computed by pppiv

uu0(~roi0) = NaN;
vv0(~roi0) = NaN;

% extend initial grid to ensure data hull covers the output grid
nr0 = length(rr0);
nc0 = length(cc0);
uu0 = [nan(1, nc0+2); nan(nr0, 1), uu0, nan(nr0, 1); nan(1, nc0+2)];
vv0 = [nan(1, nc0+2); nan(nr0, 1), vv0, nan(nr0, 1); nan(1, nc0+2)];
    
spc0 = rr0(2)-rr0(1);
rr0 = [(rr0(1)-spc0), rr0, (rr0(end)+spc0)];
cc0 = [(cc0(1)-spc0), cc0, (cc0(end)+spc0)];

% validate, replace, smooth, see [3-4])
% [uu0, vv0, sf] = pppiv(uu0, vv0);
[uu0, vv0, sf] = pppiv(uu0, vv0, 'nosmoothing');


% interpolate displacements to new sample grid
uu1 = interp2(cc0, rr0, uu0, cc1(:)', rr1(:), 'linear');
vv1 = interp2(cc0, rr0, vv0, cc1(:)', rr1(:), 'linear');

end



% verbose message subroutines --------------------------------------------------

function [] = print_sep(msg, verbose)
% Print a user-specified message and a separator line for verbose output
% messages

if verbose
    fprintf('----------\n%s\n', msg);
end

end

function print_input(ini, fin, xx, yy, npass, samplen, sampspc, umin, ...
    umax, vmin, vmax, ncbc, verbose)
% Display values (or a summary of them) for the input arguments

if verbose
    fprintf('ini: size = [%i, %i], fraction data = %.2f%%\n',...
        size(ini, 1), size(ini, 2), sum(ini(:) ~= 0)/numel(ini)*100);
    fprintf('fin: size = [%i, %i], fraction data = %.2f%%\n',...
        size(fin, 1), size(fin, 2), sum(fin(:) ~= 0)/numel(fin)*100);
    fprintf('xx: length = %i, min = %.3f, max = %.3f, delta = %.3f\n', ...
        length(xx), min(xx), max(xx), xx(2)-xx(1));
    fprintf('yy: length = %i, min = %.3f, max = %.3f, delta = %.3f\n', ...
        length(yy), min(yy), max(yy), yy(2)-yy(1));
    fprintf('npass: %i\n', npass);
    fprintf('samplen: %s\n', sprintf('%i  ', samplen));
    fprintf('sampspc: %s\n', sprintf('%i  ', sampspc));
    fprintf('umin: %s\n', sprintf('%.2f  ', umin));
    fprintf('umax: %s\n', sprintf('%.2f  ', umax));
    fprintf('vmin: %s\n', sprintf('%.2f  ', vmin));
    fprintf('vmax: %s\n', sprintf('%.2f  ', vmax));
    fprintf('ncbc: %s\n', sprintf('%i  ', ncbc));
end

end

function [] = print_pass(rr, cc, umax, umin, vmax, vmin, ncbc, verbose)
% Display information for PIV pass 
%
% Arguments:
%
% rr, cc = Vector, integer, sample grid coordinate vectors in the y/row and
%   x/column directions
%
% umax, umin, vmax, vmin = Scalar, double, displacement limits in the x-
%   and y-direction, in pixel coordinates


if verbose
    fprintf('x-dir grid: npts = %i, min = %.2f, max = %.2f\n', ...
        length(cc), min(cc), max(cc));
    fprintf('y-dir grid: npts = %i, min = %.2f, max = %.2f\n', ...
        length(rr), min(rr), max(rr));
    fprintf('u limits, pixels: min = %.2f, max = %.2f\n', ...
        umax, umin);
    fprintf('v limits, pixels, y-dir: min = %.2f, max = %.2f\n', ...
        vmax, vmin);
    fprintf('correlation-based-correction stencil = %g\n', ncbc);
end

end

function [] = plot_win_loc(data, pos)
% Plot the outline of the sample/interrogation window in pixel coordinates on
% top of the data matrix.
%
% Arguments:
% 
%   data = 2D matrix, double, data matrix
%
%   pos = Vector, double, window position vector as [left, bottom, width,
%       height] in poxel coordinates

imagesc(data)
axis equal
axis off
box off
set(gca, 'Clipping', 'off')
colormap(gray)
hold on
plot([pos(1), pos(1), pos(1)+pos(3), pos(1)+pos(3), pos(1)], ...
    [pos(2), pos(2)+pos(4), pos(2)+pos(4), pos(2), pos(2)], ...
    'Color', 'r', 'LineWidth', 2);
hold off

end

function [] = plot_xcr_plane(xwin, rpk, cpk)
% Plot the correlation plane with the location of the peak
%
% Arguments: 
%
%   xwin = 2D matrix, double, cross-correlation window
%
%   rpk, cpk = Scalar, double, location of the correlation peak in pixel
%       coordinates.

imagesc(xwin);
axis equal;
hold on
plot([1, size(xwin, 2)], [rpk, rpk], ':k');
plot([cpk, cpk], [1, size(xwin, 1)], ':k');
hold off
colorbar

end

function [] = plot_win(swin, iwin, rpk, cpk)
% Plot sample and interrogation windows, with best fit sample window footprint
%
% Arguments:
%
%   swin, iwin = 2D matrix, double, sample and interrogation windows
% 
%   rpk, cpk = Scalar, double, location of the peak in the (valid) correlation
%       plane in pixel coordinates.

subplot(1,2,1)
imagesc(swin);
axis equal
axis off
box off
colormap(gray)
title('Sample Window');

subplot(1,2,2);
imagesc(iwin);
axis equal
axis off
box off
set(gca, 'Clipping', 'off')
colormap(gray)
hold on

hwidth = (size(swin,1)-1)/2;
rpk = rpk+hwidth; % adjust indices from 'valid' extent to 'same' extent
cpk = cpk+hwidth;
rmin = rpk-hwidth;
rmax = rpk+hwidth;
cmin = cpk-hwidth;
cmax = cpk+hwidth;
plot([cmin, cmin, cmax, cmax, cmin], [rmin, rmax, rmax, rmin, rmin], ...
    'Color', 'r', 'LineWidth', 2);
hold off
title('Interrogation Window');
                
end

function [] = plot_sample_point(sdata, idata, swin, spos, iwin, ipos, xwin, rpk, cpk, verbose)
% Plot debugging information for a single sample points, pausing for user
% to advance.
%
% Arguments
% 
%   sdata, idata = 2D matrix, full data matrix for the sample (initial) and
%       interrogation (final) model states.
%
%   swin, iwin = 2D matrix, sample and interrogation windows for this sample
%       point.
%
%   spos, ipos = Vector, double, position vectors for sample and interrogation
%       windows in pixel coordinates, formatted as [left, bottom, width, height]
%
%   xwin = 2D matrix, cross-correlation plane, valid extent only
%
%   rpk, cpk = Scalar, double, sub-pixel location of the correlation plane peak
%
%   verbose = Scalar, integer, verbosity flag, must be == 2 to enable plotting

if verbose == 2 
    figure(1)
    subplot(2,1,1)
    plot_win_loc(sdata, spos)
    title('Initial State with Sample Window')
    
    figure(1)
    subplot(2,1,2)
    plot_win_loc(idata, ipos)
    title('Final State with Interrogation Window')
    
    figure(2)
    plot_xcr_plane(xwin, rpk, cpk)
    
    figure(3)
    plot_win(swin, iwin, rpk, cpk);
    
    drawnow;
    pause
end

end
