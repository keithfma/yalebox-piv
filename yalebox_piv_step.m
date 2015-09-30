function [xx, yy, uu, vv] = ...
    yalebox_piv_step(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0, verbose)
%
% function [xx, yy, uu, vv] = ...
%     yalebox_piv_step(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0. verbose)
%
% New implementation PIV analysis for Yalebox image data
%
% Arguments, input:
%
%   ini = 2D matrix, double, range 0 to 1, normalized grayscale image from
%       the start of the step to be analyzed.
%
%   fin = 2D matrix, double, range 0 to 1, normalized grayscale image from
%       the end of the step to be analyzed.
%
%   xx = Vector, double, increasing, x-direction coordinate vector, length
%       must match the columns in ini and fin.
%
%   yy = Vector, double, increasing, y-direction coordinate vector, length
%       must match the rows in ini and fin.
%
%   samplen = Scalar, integer, side length of the square sample window
%
%   sampspc = Scalar, integer, spacing between adjacent sample points in the
%       (square) sample grid
%
%   intrlen = Scalar, integer, side length of the square interrogation window
%
%   u0 = Scalar (PLANNED: or 2D matrix), integer, initial guess for x-direction
%       displacement, units are [pixels]
%
%   v0 = Scalar (PLANNED: or 2D matrix), integer, initial guess for y-direction
%       displacement, units are [pixels]
%
%   verbose = Scalar, integer, flag to enable (1) or diasable (0) verbose text
%       output messages
%
% Arguments, output:
%
%   xx, yy = Vector, double, coordinate vectors for the final output sample
%       grid, in world coordinate units
%
%   uu, vv = 2D matrix, double, computed displacement in the x- and y-directions
%       in world coordinate units
%
% References:
%
% [1] Raffel, M., Willert, C. E., Wereley, S. T., & Kompenhans, J. (2007).
%   Particle Image Velocimetry: A Practical Guide. BOOK, 1–448. 
%
% [2] Wereley, S. T. (2001). Adaptive Second-Order Accurate Particle Image
%   Velocimetry, 31
%

% parse inputs
check_input(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0, verbose);
if verbose
    print_sep('Input Arguments');
    print_input(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0);
end

% init sample grid 
[rr, cc] = sample_grid(samplen, sampspc, size(ini));
nr = length(rr);
nc = length(cc);

% init displacements
uu = u0*ones(nr, nc); % set to  initial guess
vv = v0*ones(nr, nc); 

% loop over sample grid
for jj = 1:nc
    for ii = 1:nr
   
        % get sample and interrogation windows
        samp = get_win(ini, rr(ii), cc(jj), samplen);               
        intr = get_win(fin, rr(ii), cc(jj), intrlen);
        
        % debug {
        show_samp_intr_win()
        % } debug
        
        % compute normalized cross-correlation
        
        % find integer displacement
        
    end
end

% debug {
xx = 0;
yy = 0;
% } debug 

keyboard

end

%% computational subroutines

function [] = check_input(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0, verbose)
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

validateattributes(samplen, ...
    {'numeric'}, {'scalar', 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'samplen');

validateattributes(sampspc, ...
    {'numeric'}, {'scalar', 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'sampspc');

validateattributes(intrlen, ...
    {'numeric'}, {'scalar', 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'intrlen');

validateattributes(u0, ...
    {'double'}, {'scalar'}, ...
    mfilename, 'u0');

validateattributes(v0, ...
    {'double'}, {'scalar'}, ...
    mfilename, 'u0');

validateattributes(verbose, {'numeric', 'logical'}, {'scalar', 'binary'}, ...
    mfilename, 'verbose');

end

function [rr, cc] = sample_grid(len, spc, sz)
% Create sample grid such that sample window edges fall at integer positions,
% and the grid is centered. Note that the first sample window requires 'len'
% rows and cols, and each additional window requires 'spc' rows and cols.
%
% Arguments:
%
%   len = Scalar, integer, length of sample window in pixels
%
%   spc = Scalar, integer, grid spacing in pixels
%
%   sz = Vector, length == 2, grid dimensions [rows, columns] for the
%       original input data matrices (e.g. ini and fin).
%
%   rr, cc = Vector, integer, coordinate vectors for the sample grid in the
%       y (a.k.a. row) and x (a.k.a. column) directions, in pixels

rem = mod(sz-len, spc)/2;

samp0 = floor(rem)+(len-1)/2; 
samp1 = sz-(ceil(rem)+(len-1)/2); 

rr = samp0(1):spc:samp1(1);
cc = samp0(2):spc:samp1(2);

end

function win = get_win(img, rcnt, ccnt, len)
% function win = get_win(img, rcnt, ccnt, len)
%
% Extract a sample or interrogation window from the input image, padding with
% zeros as needed.
%
% For sample windows: Sample grid centroids are chosen such that sample window
% edges lie at integer pixel coordinates. Rounding operations thus have no
% effect in this case.
%
% For interrogation windows: Window edges are not guaranteed to be integer pixel
% coordinates. Rounding is used to expand window extent outward to the nearest
% integer pixel coordinates. Note that this does not change the window centroid,
% since expansion is always symmetric.
%
% Arguments:
%
%   img = 2D matrix, double, range 0 to 1, normalized grayscale image from
%       the pair to be analyzed. Should be the initial image for the sample
%       window and the final image for the interrogation window.
%
%   rcnt, ccnt = Scalar, double, location of the window centroid 
%
%   len = Scalar, integer, length of the window

% get window limits, may lie outside the image domain
hlen = (len-1)/2;
r0 = floor(rcnt-hlen);
r1 =  ceil(rcnt+hlen); 
c0 = floor(ccnt-hlen);
c1 =  ceil(ccnt+hlen);

% get pad size, restrict window indices to valid range
pl = max(0, 1-c0);
c0 = max(1, c0);

nc = size(img, 2);
pr = max(0, c1-nc);
c1 = min(nc, c1);

pb = max(0, 1-r0);
r0 = max(1, r0);

nr = size(img, 1);
pt = max(0, r1-nr);
r1 = min(nr, r1);

% extract data and add pad
sub = img(r0:r1, c0:c1);
[snr, snc] = size(sub);
win = [zeros(pb, pl+snc+pr);
       zeros(snr, pl), sub, zeros(snr, pr);
       zeros(pt, pl+snc+pr)];  
    
end


%% verbose subroutines --------------------------------------------------

function [] = print_sep(msg)
% Print a user-specified message and a separator line for verbose output
% messages

fprintf('----------\n%s\n', msg);

end

function print_input(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0)
% Display values (or a summary of them) for the input arguments

fprintf('ini: size = [%i, %i], fraction data = %.2f%%\n',...
    size(ini, 1), size(ini, 2), sum(ini(:) ~= 0)/numel(ini)*100);

fprintf('fin: size = [%i, %i], fraction data = %.2f%%\n',...
    size(fin, 1), size(fin, 2), sum(fin(:) ~= 0)/numel(fin)*100);

fprintf('xx: length = %i, min = %.3f, max = %.3f, delta = %.3f\n', ...
    length(xx), min(xx), max(xx), xx(2)-xx(1));

fprintf('yy: length = %i, min = %.3f, max = %.3f, delta = %.3f\n', ...
    length(yy), min(yy), max(yy), yy(2)-yy(1));

fprintf('samplen: %s\n', sprintf('%i  ', samplen));

fprintf('sampspc: %s\n', sprintf('%i  ', sampspc));

fprintf('intrlen: %s\n', sprintf('%i  ', intrlen));

fprintf('u0: %i\n', u0);

fprintf('u0: %i\n', v0);

end