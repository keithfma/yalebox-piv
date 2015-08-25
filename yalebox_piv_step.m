function [] = yalebox_piv_step()
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
 
% debug { 
load('debug_input.mat', 'ini', 'fin', 'xx', 'yy');
npass = 2;
samplen = [53, 27];
yrez = [0, 0];
xrez = [50, 100];
verbose = true;
umax = [0.02, 0.01];
umin = -umax;
vmax = [0.02, 0.01];
vmin = -vmax;
validate = true;
eps0 = 1.1;
epsthresh = 2.2;
data_min_frac = 0.3;
% } debug

% parameters
nodata = 0; % value used to indicate "no data"

print_input(verbose, 'input', ini, fin, xx, yy, npass, samplen, ...
    xrez, yrez, umin, umax, vmin, vmax, validate, eps0, epsthresh, ...
    data_min_frac); 

check_input(ini, fin, xx, yy, npass, samplen, xrez, yrez, umin, umax, vmin, ...
    vmax, validate, eps0, epsthresh, data_min_frac);

[xrez, yrez] = equalize_unknown_grid_dims(xrez, yrez, size(ini,1), size(ini,2));

[umin, umax] = uv_lim_world_to_pixel(umin, umax, xx);
[vmin, vmax] = uv_lim_world_to_pixel(vmin, vmax, yy);

print_input(verbose, 'preprocessed input', ini, fin, xx, yy, npass, ...
    samplen, xrez, yrez, umin, umax, vmin, vmax, validate, eps0, epsthresh, ...
    data_min_frac); 

% % debug {
% npass = 1;
% % } debug

% initialize pixel grid and displacements for sample_grid()
rr = 1:length(yy);
cc = 1:length(xx);
uu = zeros(length(yy), length(xx));
vv = zeros(length(yy), length(xx));

% loop over PIV passes
for pp = 1:npass
    
    print_pass(verbose, pp, npass, samplen, yrez, xrez, umax, umin, vmax, vmin);
    
    [xx, yy, cc, rr, uu, vv] = sample_grid(yrez(pp), xrez(pp), xx, yy, cc, ...
                                           rr, uu, vv);
    nsamp = samplen(pp)*samplen(pp);
                                       
    % loop over sample grid    
    for ii = 1:xrez(pp)
        for jj = 1:yrez(pp)
            
            % (next: loop over correlation-based-correction samples)
            
            % get sample window
            samp = get_samp_win(ini, rr(jj), cc(ii), samplen(pp));
            
            % skip if there is too little data in the sample window
            data_frac = 1-sum(samp(:) == nodata)/nsamp;
            if data_frac < data_min_frac;
                uu(jj, ii) = nodata; % ATTN
                vv(jj, ii) = nodata; % ATTN
                continue
            end
            
            % get interrogation window
            [intr, uintr, vintr] = get_intr_win(fin, rr(jj), cc(ii), ...
                                       umin(pp)+uu(jj,ii), umax(pp)+uu(jj,ii), ...
                                       vmin(pp)+vv(jj,ii), vmax(pp)+vv(jj,ii));
            
            % compute correlation
            
            % (next: accumulate correlation-based-correction samples)
            
            % (next: end loop over correlation-based-correction samples)
            
            % find the correlation plane maximum (with subpixel accuracy)
            
            % get displacement
            
            % identify and interpolate outliers
            
        end
    end
    
end
 
end

function [] = check_input(ini, fin, xx, yy, npass, samplen, xrez, yrez, ...
                  umin, umax, vmin, vmax, validate, eps0, epsthresh, ...
                  data_min_frac)
% Check for sane input argument properties, exit with error if they do not match
% expectations.
              
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
validateattributes(validate, ...
    {'logical'}, {'scalar'}, ...
    mfilename, 'validate');
validateattributes(eps0, ...
    {'double'}, {'scalar', 'real', 'nonnegative', 'nonnan', 'finite'}, ...
    mfilename, 'eps0');
validateattributes(epsthresh, ...
    {'double'}, {'scalar', 'real', 'nonnegative', 'nonnan', 'finite'}, ...
    mfilename, 'epsthresh');
validateattributes(data_min_frac, ...
    {'double'}, {'scalar', 'nonnegative', '<=', 1, 'nonnan', 'finite'}, ...
    mfilename, 'data_min_frac');

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
% coordinate axis polarity, round to whole pixel

dxy = xy(2)-xy(1); 
uvminmax = sort([uvmin(:), uvmax(:)]/dxy, 2);
uvmin = round(uvminmax(:,1));
uvmax = round(uvminmax(:,2));

end

function [x1, y1, c1, r1, u1, v1] = sample_grid(nr1, nc1, x0, y0, c0, r0, u0, v0)
% [x1, y1, c1, r1, u1, v1] = sample_grid(nr1, nc1, x0, y0, c0, r0, u0, v0)
%
% Compute sample coordinate grid for the new pass in both pixel and world
% coordinates, and interpolate the previously computed displacements to the new
% grid. The grid is linearly spaced over the model domain, including edges.
% Output values refer to sample window center locations, which are are not in
% general integers. Note that this approach assumes the coordinate system is
% linear (can interpolate from low-res to high-res)
%
% Arguments:
%
%   nr1, nc1 = Scalars, integers, number of rows and columns in the new sample
%       grid for the new pass
%
%   x0, y0 = Vectors, x-dir (a.k.a column-dir) and y-dir (a.k.a row-dir) world
%       coordinates for the previous pass. 
%
%   c0, r0 = Vectors, x-dir (a.k.a. row-dir) and y-dir (a.k.a. row-dir) pixel
%       coordinates for the sample grid from the previous pass
%
%   u0, v0 = Matrices, x-dir and y-dir displacements computed on the sample grid
%       for the previous pass 
%
%   x1, y1, c1, r1, u1, v1 = Same as the above, but for the new pass.

% new pixel sample grid
r1 = linspace(1, max(r0), nr1);
c1 = linspace(1, max(c0), nc1);

% new world sample grid
y1 = interp1(r0, y0, r1, 'linear');
x1 = interp1(c0, x0, c1, 'linear');

% interpolate displacements from previous pass to current grid
[c1mat, r1mat] = meshgrid(c1, r1);
u1 = interp2(c0, r0, u0, c1mat, r1mat);
v1 = interp2(c0, r0, v0, c1mat, r1mat);

end

function [win] = get_samp_win(data, rpt, cpt, slen)
% Get sample window for a single sample grid point. Sample windows are squares
% with pre-defined side length, approximately centered on the sample grid point,
% zero-padded if necessary.
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

function [win, uwin, vwin] = get_intr_win(data, rpt, cpt, umin, umax, vmin, vmax)
%
% Get interrogation window for a single sample grid point. Interrogation windows
% are rectangular, with shape approximately defined by the specified minimum and
% maximum displacements. The exact displacements corresponding to each element
% in the window are returned as coordinate vectors.
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

% get window index, with range expanded to nearest whole pixel
cmin = floor(cpt+umin);
cmax = ceil(cpt+umax);
rmin = floor(rpt+vmin);
rmax = ceil(rpt+vmax);

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
win = [zeros(pt, pl+snc+pr);
       zeros(snr, pl), sub, zeros(snr, pr);
       zeros(pb, pl+snc+pr)];
   
% % debug {
% imagesc(data); caxis([0,1]); axis equal; hold on;
% plot([c0, c0, c1, c1, c0], [r0, r1, r1, r0, r0], 'LineWidth', 2, 'Color', 'k');
% hold off; drawnow; pause(0.25);
% % } debug
    
end

% verbose message subroutines --------------------------------------------------

function print_sep()
% Print a separator line for verbose output messages

fprintf('----------\n');

end

function print_input(verbose, msg, ini, fin, xx, yy, npass, samplen, ...
             xrez, yrez, umin, umax, vmin, vmax, validate, eps0, epsthresh, ...
             data_min_frac)
% Display values (or a summary of them) for the input arguments

if verbose
    print_sep;
    fprintf('%s\n', msg);
    fprintf('ini: size = [%i, %i], min = %.2f. max = %.2f, masked = %.2f%%\n',...
        size(ini, 1), size(ini, 2), min(ini(:)), max(ini(:)), ...
        sum(ini(:) == 0)/numel(ini)*100);
    fprintf('fin: size = [%i, %i], min = %.2f. max = %.2f, masked = %.2f%%\n',...
        size(fin, 1), size(fin, 2), min(fin(:)), max(fin(:)), ...
        sum(fin(:) == 0)/numel(fin)*100);
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
    fprintf('validate: %i\n', validate);
    fprintf('eps0: %.3f\n', eps0);
    fprintf('epsthresh: %.3f\n', epsthresh);
    fprintf('data_min_frac: %.3f\n', data_min_frac);
end

end

function [] = print_pass(verbose, ind, npass, samplen, yrez, xrez, umax, ...
                  umin, vmax, vmin)
% Display parameters for PIV "ind" or "npass"

if verbose
    print_sep;
    fprintf('PIV pass %i of %i\n', ind, npass);
    fprintf('samplen = %i\n', samplen(ind));
    fprintf('xrez = %i\tyrez = %i\n', xrez(ind), yrez(ind));
    fprintf('umax = %.2f\tumin = %.2f\n', umax(ind), umin(ind));
    fprintf('vmax = %.2f\tvmin = %.2f\n', vmax(ind), vmin(ind));
end

end

