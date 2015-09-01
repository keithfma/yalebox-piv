function [xx, yy, uu, vv, ss] = yalebox_piv_step()
% Re-implementation of yalebox PIV analysis routine
%
% Notes:
%
% Sample grid does not fall on integer pixels - the sample windows are
% constructed so that they are samplen pixels wide, and the grid point is as
% close to the center as possible, given that the data are discrete, the window
% center may be as much as 0.5 pixels from the grid point, this is almost
% certainly negligible, but worth noting.
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
%   samplen = Vector, length === npass, double, side length of the square sample
%       window
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
%   validate = Scalar, logical, flag to enable (true) or disable (false) vector
%       validation using the normalized median filter.
%
%   eps0 = Scalar, double, parameter to normalized median filter used for vector
%       validation, ignored if validate == false
%
%   epsthresh = " "
%
%   data_min_frac = Scalar, double, in the range [0, 1] inclusive, minimum
%       fraction of sample window that must contain data to proceed with PIV, 
%       assumes that no-data pixels are set to 0
%   
%   verbose = Scalar, integer, flag to set verbosity level, (0) no verbose
%       output, (1) enable verbose text output messages, (2) enable debugging
%       plots
%
% Arguments, output:
% 
% xx, yy, uu, vv, smoothing factors?
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

% load test case input data
load('test01_input.mat', 'ini', 'fin', 'xx', 'yy');
%load('test02_input.mat', 'ini', 'fin', 'xx', 'yy');

% single pass
npass = 1;
samplen = 30;
sampspc = 15;
verbose = 1;
umax =  0.05;
umin = -0.05;
uinit = 0;
vmax =  0.05;
vmin = -0.05;
vinit = 0;

% % dual pass
% npass = 2;
% samplen = [60, 30];
% sampspc = [30, 15];
% verbose = 1;
% umax = [ 0.02,  0.005];
% umin = [-0.02, -0.005];
% uinit = 0;
% vmax = [ 0.02,  0.005];
% vmin = [-0.02, -0.005];
% vinit = 0;

% % tri pass
% npass = 3;
% samplen = [60, 30, 15];
% sampspc = [30, 15, 7];
% verbose = 1;
% umax = [ 0.02,  0.005,  0.0025];
% umin = [-0.02, -0.005, -0.0025];
% uinit = 0;
% vmax = [ 0.02,  0.005,  0.0025];
% vmin = [-0.02, -0.005, -0.0025];
% vinit = 0;

% } debug

print_sep('INPUT ARGUMENTS', verbose);
print_input(ini, fin, xx, yy, npass, samplen, sampspc, umin, umax, ...
        uinit, vmin, vmax, vinit, verbose);

check_input(ini, fin, xx, yy, npass, samplen, sampspc, umin, umax, ...
    vmin, vmax);

[umin, umax, uinit] = uv_input_world_to_pixel(umin, umax, uinit, xx);
[vmin, vmax, vinit] = uv_input_world_to_pixel(vmin, vmax, vinit, yy);

[nr0, nc0] = size(ini);
[rr, cc] = sample_grid(sampspc(1), nr0, nc0);
uu = uinit*ones(length(rr), length(cc));
vv = vinit*ones(length(rr), length(cc));

% loop over PIV passes
ss = nan(npass, 1);
for pp = 1:npass
    
    print_sep(sprintf('PIV pass %i of %i', pp, npass), verbose);
    print_pass(rr, cc, umax(pp), umin(pp), vmax(pp), vmin(pp), verbose);
    
    roi = true(length(rr), length(cc));
                                       
    % loop over sample grid    
    for ii = 1:length(cc)
        for jj = 1:length(rr)
            
            % (next: loop over correlation-based-correction samples)
            
            % skip if window center lies outside the roi at start or finish
            rchk = round(rr(jj));
            cchk = round(cc(ii));            
            if ini(rchk, cchk) == 0 || fin(rchk,cchk) == 0
                roi(jj, ii) = false;
                continue
            end
            
            % get sample window
            [samp, srmin, srmax, scmin, scmax]  = ...
                get_samp_win(ini, rr(jj), cc(ii), samplen(pp));
            
            % get interrogation window
            [intr, uintr, vintr, irmin, irmax, icmin, icmax] = ...
                get_intr_win(fin, rr(jj), cc(ii), ...
                    umin(pp)+uu(jj,ii), umax(pp)+uu(jj,ii), ...
                    vmin(pp)+vv(jj,ii), vmax(pp)+vv(jj,ii), samplen(pp));
                            
            % skip if interrogation window is empty
            if max(intr(:)) == 0
                roi(jj, ii) = false;
                continue
            end
                
            % compute correlation, trimming to valid range (see help)
            xcr = get_cross_corr(samp, intr);
                                    
            % (next: accumulate correlation-based-correction samples)
            
            % (next: end loop over correlation-based-correction samples)           
            
            % find the correlation plane maximum with subpixel accuracy
            [rpeak, cpeak] = find_peak(xcr);
            
            % get displacement in pixel coordinates
            vv(jj, ii) = interp1(1:size(xcr, 1), vintr, rpeak);    
            uu(jj, ii) = interp1(1:size(xcr, 2), uintr, cpeak);
            
            plot_sample_point(ini, fin, samp, srmin, srmax, scmin, ...
                scmax, intr, irmin, irmax, icmin, icmax, xcr, rpeak, ...
                cpeak, verbose)
                    
        end
    end
        
    
    % post-process and prep for next pass
    pp_next = min(npass, pp+1); % last pass uses same grid
    [rr_next, cc_next] = sample_grid(sampspc(pp_next), nr0, nc0);            
%     [uu, vv] = post_process(uu, vv, roi, rr, cc, rr_next, cc_next);
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
 
% % debug {
% keyboard
% % } debug

end

function check_input(ini, fin, xx, yy, npass, samplen, sampspc, umin, ...
             umax, vmin, vmax)
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

end

function [uvmin, uvmax, uvinit] = ...
	uv_input_world_to_pixel(uvmin, uvmax, uvinit, xy)
% Convert displacement limits from world to pixel coordinates, one
% direction at a time. Sort to preserve the correct min/max regardless of
% the world coordinate axis polarity, round to whole pixel.
%
% Arguments:
%
%   uvmin, uvmax = Scalar, double, minimum and maximum displacement in
%       world coordinates for either x- or y-direction.
%
%   uvinit = Scalar, double, initial guess for displacement in world
%       coordinates for either x- or y- direction.
%
%   xy = Vector, double, world coordinate vector for either x- or
%       y-direction

dxy = xy(2)-xy(1); 

uvminmax = sort([uvmin(:), uvmax(:)]/dxy, 2);
uvmin = round(uvminmax(:,1));
uvmax = round(uvminmax(:,2));

uvinit = uvinit/dxy;

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

function [win, rmin, rmax, cmin, cmax] = get_samp_win(data, rpt, cpt, slen)
% Get sample window for a single sample grid point. Sample windows are
% squares with pre-defined side length, approximately centered on the
% sample grid point, zero-padded if necessary.
%
% Arguments:
%
%   data = 2D Matrix, initial model state data
%
%   rpt, cpt = Scalar, double, row-, col-position of the sample point (window
%       center)
%
%   slen = Scalar, integer, side length of square sample window
%
%   win = 2D matrix, double, sample window containing a subset of the input data
%       and perhaps some zero padding
%
%   rmin, rmax, cmin, cmax = Scalar, integer, window limits in pixel
%       coordinates

% get window index range shifted to nearest whole pixel
rmin = rpt-(slen-1)/2;
radj = round(rmin)-rmin;
rmin = rmin+radj;
rmax = rpt+(slen-1)/2+radj;

cmin = cpt-(slen-1)/2;
cadj = round(cmin)-cmin;
cmin = cmin +cadj;
cmax = cpt+(slen-1)/2+cadj;

% extract subset, including pad if needed
win = get_padded_subset(data, rmin, rmax, cmin, cmax);

end

function [win, uwin, vwin, rmin, rmax, cmin, cmax] = ...
    get_intr_win(data, rpt, cpt, umin, umax, vmin, vmax, slen)
%
% Get interrogation window for a single sample grid point. Interrogation windows
% are rectangular, with shape approximately defined by the specified minimum and
% maximum displacements. The exact displacements corresponding to each element
% in the window are returned as coordinate vectors. A border the size of
% half the sample window side length is added on all sides to ensure that
% all the tested displacements lie within the "valid" range of the cross
% correlation output.
%
% Arguments:
%
%   data = 2D Matrix, final model state data
%
%   rpt, cpt = Scalar, double, row-, col-position of the sample point (window
%       center)
%
%   umin, umax = Scalar, double, minimum, maximum x-direction displacements in
%       pixel coordinates
%
%   vmin, vmax = Scalar, double, minimum, maximum y-direction displacements in
%       pixel coordinates
%
%   win = 2D matrix, double, interrogation window containing a subset of the
%       input data and perhaps some zero padding
%
%   uwin = Vector, length == size(win,2), coordinate vector for the
%       interrogation window giving x-direction displacements in pixel
%       coordinates.
%
%   vwin = Vector, length == size(win,1), coordinate vector for the
%       interrogation window giving y-direction displacements in pixel
%       coordinates.
%
%   slen = Scalar, integer, side length of square sample window
%
%   rmin, rmax, cmin, cmax = Scalar, integer, window limits in pixel
%       coordinates

% get size of boundary to include (sample window half-width)
bnd = (slen-1)/2;

% get window index, with range expanded to nearest whole pixel
cmin = floor(cpt+umin-bnd);
cmax = ceil(cpt+umax+bnd);
rmin = floor(rpt+vmin-bnd);
rmax = ceil(rpt+vmax+bnd);

% get displacements from sample center to rows and cols of interrogation window
uwin = (cmin:cmax)-cpt;
vwin = (rmin:rmax)-rpt;

% extract subset, including pad if needed
win = get_padded_subset(data, rmin, rmax, cmin, cmax);

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
% Compute normalized cross correlation, and crop to the extent of larger
% matrix (bb). Allows for non-square aa, although that is not needed at
% this time.
%
% Arguments:
%   aa = 2D matrix, double, smaller 'template' matrix, as used here, this
%       is the sample window
%
%   bb = 2D matrix, double, larger matrix, as used here, this is the
%       interrogation windo

fullxcorr = normxcorr2(aa, bb);

% compute pad size in both dimensions
aaSize = size(aa);
npre = floor(aaSize/2); % pre-pad
npost = aaSize-npre-1; % post-pad

xcorr = fullxcorr( (1+npre(1)):(end-npost(1)), (1+npre(2)):(end-npost(2)) );
            
end

function [rpk, cpk] = find_peak(zz)
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

[rpk, cpk] = find(zz == max(zz(:)));

% peak is at the edge of the matrix, cannot compute subpixel location
if rpk == 1 || rpk == size(zz, 1) || cpk == 1 || cpk == size(zz,2)
    rpk = -1;
    cpk = -1;
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
    rpk = -1;
    cpk = -1;
    
end

end

function [uu1, vv1] = post_process(uu0, vv0, roi0, rr0, cc0, rr1, cc1)
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
[uu0, vv0, sf] = pppiv(uu0, vv0);

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
    umax, uinit, vmin, vmax, vinit, verbose)
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
    fprintf('uinit: %s\n', sprintf('%.2f  ', uinit));
    fprintf('vmin: %s\n', sprintf('%.2f  ', vmin));
    fprintf('vmax: %s\n', sprintf('%.2f  ', vmax));
    fprintf('vinit: %s\n', sprintf('%.2f  ', vinit));
end

end

function [] = print_pass(rr, cc, umax, umin, vmax, vmin, verbose)
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
end

end

function [] = plot_win_loc(data, rmin, rmax, cmin, cmax)
% Plot the data matrix with the sample window outlined in red

imagesc(data)
axis equal
axis off
box off
set(gca, 'Clipping', 'off')
colormap(gray)
hold on
plot([cmin, cmin, cmax, cmax, cmin], [rmin, rmax, rmax, rmin, rmin], ...
    'Color', 'r', 'LineWidth', 2);
hold off

end

function [] = plot_xcr_plane(xcr, rpeak, cpeak)
% Plot the correlation plane with the location of the peak

imagesc(xcr);
axis equal;
hold on
plot([1, size(xcr, 2)], [rpeak, rpeak], ':k');
plot([cpeak, cpeak], [1, size(xcr, 1)], ':k');
hold off
colorbar

end

function [] = plot_win(samp, intr, rpeak, cpeak, slen)
% Plot sample and interrogation windows, with best fit sample window footprint

subplot(1,2,1)
imagesc(samp);
axis equal
axis off
box off
colormap(gray)
title('Sample Window');

subplot(1,2,2);
imagesc(intr);
axis equal
axis off
box off
set(gca, 'Clipping', 'off')
colormap(gray)
hold on
rmin = rpeak-(slen-1)/2;
rmax = rpeak+(slen-1)/2;
cmin = cpeak-(slen-1)/2;
cmax = cpeak+(slen-1)/2;
plot([cmin, cmin, cmax, cmax, cmin], [rmin, rmax, rmax, rmin, rmin], ...
    'Color', 'r', 'LineWidth', 2);
hold off
title('Interrogation Window');
                
end

function [] = plot_sample_point(ini, fin, samp, srmin, srmax, scmin, ...
                  scmax, intr, irmin, irmax, icmin, icmax, xcr, rpeak, ...
                  cpeak, verbose)
% Plot debugging information for a single sample points, pausing for user
% to advance.
%
% Arguments
% 
% 

if verbose == 2 
    figure(1)
    subplot(2,1,1)
    plot_win_loc(ini, srmin, srmax, scmin, scmax)
    title('Initial State with Sample Window')
    
    figure(1)
    subplot(2,1,2)
    plot_win_loc(fin, irmin, irmax, icmin, icmax)
    title('Final State with Interrogation Window')
    
    figure(2)
    plot_xcr_plane(xcr, rpeak, cpeak)
    
    figure(3)
    plot_win(samp, intr, rpeak, cpeak, size(samp,1));
    
    drawnow;
    pause
end

end
