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
%   Particle Image Velocimetry: A Practical Guide. BOOK, 1???448. 
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
        [samp, samp_pos] = get_win(ini, rr(ii), cc(jj), samplen); 
        [intr, intr_pos] = get_win(fin, rr(ii), cc(jj), intrlen); 
        
        % % debug {
        % figure(1)
        % show_win(ini, fin, rr(ii), cc(jj), samp, samp_pos, intr, intr_pos);
        % % } debug
        
        % compute normalized cross-correlation
        xcr = normxcorr2(samp, intr);
        
        % % find correlation plane max, integer pixel precision
        % [rpeak, cpeak] = find(xcr == max(xcr(:)));
        
        % find correlation plane max, subpixel precision
        [rpeak, cpeak] = get_peak_centroid(xcr);
        
        % find displacement from position of the correlation max
        %   - account for padding (-samplen)
        %   - account for relative position of interogation and sample
        %     windows (e,g, for columns: -(samp_pos(1)-intr_pos(1))         
        uu(ii, jj) = cpeak-samplen-(samp_pos(1)-intr_pos(1));
        vv(ii, jj) = rpeak-samplen-(samp_pos(2)-intr_pos(2));
        
    end
end

% convert displacements to world coordinates (assumes constant grid spacing)
uu = uu.*(xx(2)-xx(1));
vv = vv.*(yy(2)-yy(1));

% interpolate world coordinates for displacement vectors
xx = interp1(1:size(ini,2), xx, cc);
yy = interp1(1:size(ini,1), yy, rr);

% debug {
keyboard
% } debug

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

function [rvec, cvec] = sample_grid(len, spc, sz)
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
%   rvec, cvec = Vector, integer, coordinate vectors for the sample grid in the
%       y (a.k.a. row) and x (a.k.a. column) directions, in pixels

rem = mod(sz-len, spc)/2;

samp0 = floor(rem)+(len-1)/2; 
samp1 = sz-(ceil(rem)+(len-1)/2); 

rvec = samp0(1):spc:samp1(1);
cvec = samp0(2):spc:samp1(2);

end

function [win, pos] = get_win(img, rcnt, ccnt, len)
% function [win, pos] = get_win(img, rcnt, ccnt, len)
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
%
%   win = 2D matrix, double, subset of img, possibly with zero padding
%
%   pos = Vector, length == 4, position vector for 'win' in the format
%       [left, bottom, width, height] in pixel coordinates
% %

% get window limits, may lie outside the image domain
hlen = (len-1)/2;
r0 = floor(rcnt-hlen);
r1 =  ceil(rcnt+hlen); 
c0 = floor(ccnt-hlen);
c1 =  ceil(ccnt+hlen);

% generate position vector for output
pos = [c0, r0, c1-c0+1, r1-r0+1];

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

function [rpk, cpk] = get_peak_centroid(xcor)
%
% Find the location of the correlation plane maximum with subpixel
% precision using a 9-point centroid calculation.
%
% Arguments:
%
%   xcor
%
%   rpk, cpk
% %

% peak location with integer precision
[rpk_int, cpk_int] = find(xcor == max(xcor(:)));

% extract subscripts and valuesfor 9-point neighborhood, trim to xcor edges 
rpk_nbr = max(1, rpk_int-1):min(size(xcor,1), rpk_int+1); 
cpk_nbr = max(1, cpk_int-1):min(size(xcor,2), cpk_int+1);
xcor_nbr = xcor(rpk_nbr, cpk_nbr);
[cpk_nbr, rpk_nbr] = meshgrid(cpk_nbr, rpk_nbr);

% compute centroid
rpk = sum(rpk_nbr(:).*xcor_nbr(:))/sum(xcor_nbr(:)); 
cpk = sum(cpk_nbr(:).*xcor_nbr(:))/sum(xcor_nbr(:));

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

function show_win(img0, img1, rcnt, ccnt, swin, spos, iwin, ipos, sec) %#ok!
% function show_win(img0, img1, rcnt, ccnt, swin, spos, iwin, ipos, sec)
%
% Display the sample and interrogation windows, as well as their position
%
% Arguments:
%
%   img0, img1 = 2D matrix, double, initial and final images
%
%   swin, iwin = 2D matrix, double, sample and interrogation windows
%
%   spos, ipos = Vector, length == 4, integer, position vectors for sample
%       and interrogation windows, formatted as [left, bottom, width,
%       height] in pixel coordinates
%
%   sec = Scalar, double, duration to pause between samples, default is to
%       wait for user button press
% %

% set defaults
if nargin < 9
    sec = [];
end

% init figure
clim = [min(img0(:)), max(img0(:))];
set(gcf, 'units', 'normalized', 'position', [0.05, 0.05, 0.9, 0.9]);

% plot initial image with window positions superimposed
subplot(2, 2, 1);
imagesc(img0);
set(gca, 'YDir', 'normal');
caxis(clim);
hold on
plot(ccnt, rcnt, 'Color', 'k', 'Marker', '*')
plot([spos(1), spos(1)+spos(3)-1, spos(1)+spos(3)-1, spos(1)          , spos(1)], ...
     [spos(2), spos(2)          , spos(2)+spos(4)-1, spos(2)+spos(4)-1, spos(2)], ...
     'Color', 'k', 'LineWidth', 2, 'LineStyle', '-');
plot([ipos(1), ipos(1)+ipos(3)-1, ipos(1)+ipos(3)-1, ipos(1)          , ipos(1)], ...
     [ipos(2), ipos(2)          , ipos(2)+ipos(4)-1, ipos(2)+ipos(4)-1, ipos(2)], ...
     'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
title('initial image');
legend({'center', 'sample', 'interrogation'}, 'Location', 'NorthEast');
hold off
axis equal
axis tight

% plot final image with window positions superimposed
subplot(2, 2, 2);
imagesc(img1);
set(gca, 'YDir', 'normal');
caxis(clim);
hold on
plot(ccnt, rcnt, 'Color', 'k', 'Marker', '*')
plot([spos(1), spos(1)+spos(3)-1, spos(1)+spos(3)-1, spos(1)          , spos(1)], ...
     [spos(2), spos(2)          , spos(2)+spos(4)-1, spos(2)+spos(4)-1, spos(2)], ...
     'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([ipos(1), ipos(1)+ipos(3)-1, ipos(1)+ipos(3)-1, ipos(1)          , ipos(1)], ...
     [ipos(2), ipos(2)          , ipos(2)+ipos(4)-1, ipos(2)+ipos(4)-1, ipos(2)], ...
     'Color', 'k', 'LineWidth', 2, 'LineStyle', '-');
title('final image');
legend({'center', 'sample', 'interrogation'}, 'Location', 'NorthEast');
hold off
axis equal
axis tight


% plot sample window
subplot(2, 2, 3);
imagesc(swin);
set(gca, 'YDir', 'normal');
caxis(clim);
title('sample window')
axis equal
axis tight
grid on

% plot interrogation window
subplot(2, 2, 4);
imagesc(iwin);
set(gca, 'YDir', 'normal');
caxis(clim);
title('interrogation window')
axis equal
axis tight
grid on

% wait for user, or pause for predetermined duration
if isempty(sec)
    pause
else 
    pause(sec)
end
 
end