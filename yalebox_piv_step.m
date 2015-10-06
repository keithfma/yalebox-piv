function [xx, yy, uu, vv] = ...
    yalebox_piv_step(ini, fin, xx, yy, samplen, sampspc, intrlen, ...
                     npass, valid_max, valid_eps, verbose)                 
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
%   samplen = Vector, length == number of grid resolutions, integer, side
%       length of the square sample window
%
%   sampspc = Vector, length == number of grid resolutions, integer,
%       spacing between adjacent sample points in the (square) sample grid
%
%   intrlen = Vector, length == number of grid resolutions, integer, side
%       length of the square interrogation window
%
%   npass = Vector, length == number of grid resolutions, integer, number of
%       image deformation passes
%
%   valid_max = Scalar, double, maximum value for the normalized residual
%       in the vector validation function, above which a vector is flagged
%       as invalid. Ref [3] reccomends a value of 2.
%
%   epsilon = Scalar, double, minumum value of the normalization factor in
%       the vector validation function. Ref [3] reccomends a value of 0.1.
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
%   Particle Image Velocimetry: A Practical Guide. BOOK. 
%
% [2] Wereley, S. T. (2001). Adaptive Second-Order Accurate Particle Image
%   Velocimetry, 31
%
% [3] Westerweel, J., & Scarano, F. (2005). Universal outlier detection for PIV
%   data. Experiments in Fluids, 39(6), 1096???1100. doi:10.1007/s00348-005-0016-6

warning off all

% parse inputs
check_input(ini, fin, xx, yy, samplen, sampspc, intrlen, npass, ...
    valid_max, valid_eps, verbose);

% report input arguments
if verbose
    print_sep('Input Arguments');
    print_input(ini, fin, xx, yy, samplen, sampspc, intrlen, npass, ...
        valid_max, valid_eps);
end

% init coordinates for full image
[cc0, rr0] = meshgrid(1:size(ini,2), 1:size(ini, 1));

% init sample grid 
[rr, cc] = sample_grid(samplen(1), sampspc(1), size(ini));
nr = length(rr);
nc = length(cc);

% init displacements as zero
uu = zeros(nr, nc); 
vv = zeros(nr, nc); 

% loop over grid resolutions
ngrid = length(samplen);
for gg = 1:ngrid
    
    % report grid refinement step
    if verbose
        print_sep(sprintf('Grid refinement step %i of %i', gg, ngrid));
    end
    
    % copy old sample grid
    rr_old = rr;
    cc_old = cc;
    
    % get new sample grid    
    [rr, cc] = sample_grid(samplen(gg), sampspc(gg), size(ini));
    nr = length(rr);
    nc = length(cc);
    
    % interpolate/extrapolate displacements to new sample grid    
    [cc_g_old, rr_g_old] = meshgrid(cc_old, rr_old);
    [cc_g,     rr_g]     = meshgrid(cc    , rr    );    
    interpolant = griddedInterpolant(rr_g_old, cc_g_old, uu, 'spline', 'spline');
    uu = interpolant(rr_g, cc_g);
    interpolant.Values = vv;
    vv = interpolant(rr_g, cc_g);
    
    % loop over image deformation passes
    for pp = 1:npass(gg)
        
        % report image deformation step
        if verbose
            print_sep(sprintf('Image deformation pass %i of %i', pp, npass(gg)));
        end
        
        % interpolate/extrapolate displacement vectors to full image resolution
        interpolant = griddedInterpolant(repmat(rr(:), 1, nc), ...
            repmat(cc(:)', nr, 1), uu, 'spline', 'spline');
        uu0 = interpolant(rr0, cc0);
        interpolant.Values = vv;
        vv0 = interpolant(rr0, cc0);
        
        % deform images (does nothing if uu0 and vv0 are 0)
        defm_ini = imwarp(ini, -cat(3, uu0, vv0)/2, 'cubic',...
            'FillValues', 0); %, 'SmoothEdges', false);
        defm_fin = imwarp(fin,  cat(3, uu0, vv0)/2, 'cubic', ...
            'FillValues', 0); %, 'SmoothEdges', false);
               
        % set mask to true
        mask = true(nr, nc);
        
        % loop over sample grid
        for jj = 1:nc
            for ii = 1:nr
                
                % get sample and (offset) interrogation windows
                [samp, samp_pos] = get_win(defm_ini, rr(ii), cc(jj), samplen(gg));
                [intr, intr_pos] = get_win(defm_fin, rr(ii), cc(jj), intrlen(gg));
                
                % skip if:
                %   - sample window is <75% full
                %   - interrogation window is empty
                if sum(samp(:) == 0) > 0.25*samplen(gg)^2 || all(intr(:) == 0)
                    
                    uu(ii, jj) = NaN;
                    vv(ii, jj) = NaN;
                    mask(ii, jj) = false;
                    continue
                end
                
                % compute normalized cross-correlation
                xcr = normxcorr2(samp, intr);
                
                % find correlation plane max, subpixel precision
                [rpeak, cpeak, stat] = get_peak_gauss2d(xcr);
                if stat == false
                    uu(ii, jj) = NaN;
                    vv(ii, jj) = NaN;
                    continue
                end
                
                % find displacement from position of the correlation max
                %   - account for padding (-samplen(gg))
                %   - account for relative position of interogation and sample
                %     windows (e,g, for columns: -(samp_pos(1)-intr_pos(1))
                delta_uu = cpeak-samplen(gg)-(samp_pos(1)-intr_pos(1));
                delta_vv = rpeak-samplen(gg)-(samp_pos(2)-intr_pos(2));
                
                uu(ii, jj) = uu(ii, jj)+delta_uu;
                vv(ii, jj) = vv(ii, jj)+delta_vv;
                
                % % debug {
                % figure(1)
                % show_win(defm_ini, defm_fin, rr(ii), cc(jj), samp, samp_pos, intr, intr_pos);
                % figure(2)
                % show_xcor(xcr, rpeak, cpeak);
                % pause
                % % } debug
                
            end % ii
        end % jj
        
        % find and drop invalid displacement vectors
        drop = validate_normalized_median(uu, vv, valid_max, valid_eps);        
        uu(drop) = NaN;
        vv(drop) = NaN;
        
        % validate, smooth, and interpolate (DCT-PLS)
        [uu, vv] = pppiv(uu, vv);
        
        % % debug {
        % keyboard
        % % } debug
        
    end % pp
    
end % gg

% re-apply mask
uu(~mask) = NaN;
vv(~mask) = NaN;

% convert displacements to world coordinates (assumes constant grid spacing)
uu = uu.*(xx(2)-xx(1));
vv = vv.*(yy(2)-yy(1));

% interpolate world coordinates for displacement vectors
xx = interp1(1:size(ini,2), xx, cc);
yy = interp1(1:size(ini,1), yy, rr);

% % debug {
% keyboard
% % } debug

end



%% computational subroutines

function [] = check_input(ini, fin, xx, yy, samplen, sampspc, intrlen, ...
                  npass, valid_max, valid_eps, verbose)
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
    {'numeric'}, {'vector', 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'samplen');

ngrid = numel(samplen);

validateattributes(sampspc, ...
    {'numeric'}, {'vector', 'numel', ngrid, 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'sampspc');

validateattributes(intrlen, ...
    {'numeric'}, {'vector', 'numel', ngrid, 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'intrlen');

validateattributes(npass, ...
    {'numeric'}, {'vector', 'numel', ngrid, 'integer', 'positive'}, ...
    mfilename, 'npass');

validateattributes(valid_max, ...
    {'double'}, {'scalar', 'positive'}, ...
    mfilename, 'valid_max');

validateattributes(valid_eps, ...
    {'double'}, {'scalar', 'positive'}, ...
    mfilename, 'valid_eps');

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

function [rpk, cpk, stat] = get_peak_gauss2d(zz)
% Find the position of the peak in matrix zz with subpixel accuracy. Peak
% location is determined from an explicit solution of two-dimensional
% Gaussian regression (REF). If the peak cannot be fit at subpixel
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

function invalid = validate_normalized_median(uu, vv, max_norm_res, epsilon)
%
% Validate the displacement vector field using a normalized median test. See
% reference [3] for details. 
%
% EXPERIMENT: increase the window size from 3x3 to 5x5
%
% Arguments:
%
%   uu, vv = 2D matrix, double, displacement vector components. NaNs are
%       treated as missing values to allow for roi masking.
%
%   max_norm_res = Scalar, double, maximum value foe the normalized residual,
%       above which a vector is flagged as invalid. Refernence [3] reccomends a
%       value of 2.
%
%   epsilon = Scalar, double, minumum value of the normalization factor.
%       Reference [3] reccomends a value of 0.1.
%
%   invalid = 2D matrix, logical, flags identifying all invalid vectors as 1
% %

% set defaults
if nargin < 3; max_norm_res = 2; end
if nargin < 4; epsilon = 0.1; end 

% init
[nr, nc] = size(uu);
invalid = false(nr, nc);
roffset = [ 2,  2,  2,  2,  2, ...
            1,  1,  1,  1,  1, ...
            0,  0,      0,  0, ...
           -1, -1, -1, -1, -1, ...
           -2, -2, -2, -2, -2];
coffset = [-2, -1,  0,  1,  2, ...
           -2, -1,  0,  1,  2, ...
           -2, -1,      1,  2, ...
           -2, -1,  0,  1,  2, ...
           -2, -1,  0,  1,  2];


% loop over all displacement vectors
for ii = 1:nr
    for jj = 1:nc
        
        % get linear indices of 8 (or less) neighbors
        rnbr = max(1, min(nr, ii+roffset));
        cnbr = max(1, min(nc, jj+coffset));
        knbr = rnbr+(cnbr-1)*nr;
        
        % extract displacements for center and neighbors
        u0 = uu(ii, jj);
        v0 = vv(ii, jj);
        unbr = uu(knbr);
        vnbr = vv(knbr);
        
        % compute neighbor median, residual, and median residual 
        med_unbr = nanmedian(unbr);
        res_unbr = abs(unbr-med_unbr);
        med_res_unbr = nanmedian(res_unbr);
        
        med_vnbr = nanmedian(vnbr);
        res_vnbr = abs(vnbr-med_vnbr);
        med_res_vnbr = nanmedian(res_vnbr);
        
        % compute center normalized residual
        norm_res_u0 = abs(u0-med_unbr)/(med_res_unbr+epsilon);
        norm_res_v0 = abs(v0-med_vnbr)/(med_res_vnbr+epsilon);
        
        % combine vector components (max or sum)
        norm_res = max(norm_res_u0, norm_res_v0);
        
        % classify as valid or invalid
        invalid(ii, jj) = norm_res > max_norm_res;
        
    end
end

end

%% verbose subroutines --------------------------------------------------

function [] = print_sep(msg)
% Print a user-specified message and a separator line for verbose output
% messages

fprintf('----------\n%s\n', msg);

end

function print_input(ini, fin, xx, yy, samplen, sampspc, intrlen, ...
                     npass, valid_max, valid_eps)
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

fprintf('npass: %i\n', npass);

fprintf('valid_max: %f\n', valid_max);

fprintf('valid_eps: %f\n', valid_eps);

end

function show_win(img0, img1, rcnt, ccnt, swin, spos, iwin, ipos) %#ok!
% function show_win(img0, img1, rcnt, ccnt, swin, spos, iwin, ipos)
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
% %


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
 
end

function [] = show_xcor(xcor, rpk, cpk) %#ok!
% plot correlation plane with the position of the peak

imagesc(xcor);
set(gca, 'YDir', 'normal');
caxis([-1 1]);
colorbar
hold on
plot(cpk, rpk, 'Color', 'k', 'Marker', '*')
title('cross-correlation');
hold off
axis equal
axis tight
end