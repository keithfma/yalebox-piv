function [xx, yy, uu, vv] = ...
    yalebox_piv(ini_ti, fin_tf, ini_roi_ti, fin_roi_tf, xx, yy, samplen, sampspc, ...
        intrlen, npass, valid_max, valid_eps, verbose)                 
% New implementation PIV analysis for Yalebox image data
%
% Arguments, input:
%
%   ini_ti, fin_tf = 2D matrix, double, range 0 to 1, normalized grayscale image from
%       the start and end of the step to be analyzed.
%
%   ini_roi_ti, fin_roi_tf = 2D matrix, logical, mask indicating pixels where there is
%       sand (1) and where there is only background (0) that should be ignored.
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

% local parameters
min_frac_data = 0.5;
min_frac_overlap = min_frac_data/2;
tension = 0.95;
roi_epsilon = 1e-2; % numerical precision for roi deformation
        
% parse inputs
check_input(ini_ti, fin_tf, ini_roi_ti, fin_roi_tf, xx, yy, samplen, sampspc, intrlen, ...
    npass, valid_max, valid_eps, verbose);

% expand grid definition vectors to reflect the number of passes
[samplen, intrlen, sampspc] = expand_grid_def(samplen, intrlen, sampspc, npass);

% init coordinate grids
[rvec, cvec] = yalebox_piv_sample_grid(samplen(1), sampspc(1), size(ini_ti));
[cc_p_tm, rr_p_tm] = meshgrid(cvec, rvec);
[cc_i, rr_i] = meshgrid(1:size(ini_ti, 2), 1:size(ini_ti ,1));

% init displacement matrices
sz = size(cc_p_tm);
uu_p_tm = zeros(sz); 
vv_p_tm = zeros(sz); 
uu_c_tm = zeros(sz); 
vv_c_tm = zeros(sz); 

sz = size(ini_ti);
uu_i_ti = zeros(sz);
vv_i_ti = zeros(sz);
uu_i_tf = zeros(sz);
vv_i_tf = zeros(sz);

% init deformed images
ini_tm = ini_ti;
fin_tm = fin_tf;

% % multipass loop
np = length(samplen)-1; 
for pp = 1:np
    
    % reset per-pass variables
    sz = size(cc_p_tm);
    cc_c_tm = zeros(sz);
    rr_c_tm = zeros(sz);
    uu_c_tm = nan(sz);
    vv_c_tm = nan(sz);
    min_overlap = min_frac_overlap*samplen(pp)*samplen(pp);    
    
    % get corrector displacements on the predictor grid
    for kk = 1:numel(uu_c_tm)
            
        % get sample window and its centroid location
        [samp, samp_pos, frac_data, r_samp_cntr, c_samp_cntr] = ...
            yalebox_piv_window(ini_tm, rr_p_tm(kk), cc_p_tm(kk), samplen(pp));
        
        % get interrogation window
        [intr, intr_pos] = ...
            yalebox_piv_window(fin_tm, rr_p_tm(kk), cc_p_tm(kk), intrlen(pp));
        
        % skip if sample window is too empty
        if frac_data < min_frac_data
            uu_c_tm(kk) = NaN;
            vv_c_tm(kk) = NaN;
            continue
        end
        
        % compute masked normalized cross-correlation
        [xcr, overlap] = normxcorr2_masked(intr, samp, intr~=0, samp~=0);
        xcr = xcr.*double(overlap>min_overlap);
        
        % find correlation plane max, subpixel precision (failed pixels -> NaN)
        [rpeak, cpeak] = yalebox_piv_peak_gauss2d(xcr);
        
        % convert position of the correlation max to displacement
        %   - account for padding in cross-correlation (-samplen(gg))
        %   - account for relative position of interogation and sample
        %     windows (e,g, for columns: -(samp_pos(1)-intr_pos(1))
        uu_c_tm(kk) = cpeak-samplen(pp)-(samp_pos(1)-intr_pos(1));
        vv_c_tm(kk) = rpeak-samplen(pp)-(samp_pos(2)-intr_pos(2));
        
        % compute location of this observation at midpoint time
        cc_c_tm(kk) = c_samp_cntr + 0.5*uu_c_tm(kk);
        rr_c_tm(kk) = r_samp_cntr + 0.5*vv_c_tm(kk);
        
    end 
    
    % interpolate corrector points to predictor grid, update predictor 
    from = ~isnan(uu_c_tm); 
    to = true(size(uu_c_tm));
    uu_c_tm(to) = spline2d(cc_p_tm(to), rr_p_tm(to), cc_c_tm(from), rr_c_tm(from), uu_c_tm(from), tension);
    vv_c_tm(to) = spline2d(cc_p_tm(to), rr_p_tm(to), cc_c_tm(from), rr_c_tm(from), vv_c_tm(from), tension);    
    
    % update predictor
    uu_p_tm = uu_p_tm + uu_c_tm;
    vv_p_tm = vv_p_tm + vv_c_tm;
    
    % validate predictor vectors (invalid -> NaN)
    [uu_p_tm, vv_p_tm] = yalebox_piv_valid_nmed(uu_p_tm, vv_p_tm, 8, valid_max, valid_eps);    
    
    % interpolate/extrapolate invalid displacement vectors
    from = ~isnan(uu_p_tm);
    to = ~from;
    uu_p_tm(to) = spline2d(cc_p_tm(to), rr_p_tm(to), cc_p_tm(from), rr_p_tm(from), uu_p_tm(from), tension);
    vv_p_tm(to) = spline2d(cc_p_tm(to), rr_p_tm(to), cc_p_tm(from), rr_p_tm(from), vv_p_tm(from), tension);
    
    % smooth predictors, 3x3 kernel smoother, with NaNs to delete points affected
    % ...by the boundary.
    uu_p_tm = padarray(uu_p_tm, [1 1], NaN, 'both');    
    vv_p_tm = padarray(vv_p_tm, [1 1], NaN, 'both');        
    kernel = fspecial('average', 3);    
    uu_p_tm = conv2(uu_p_tm, kernel, 'same');
    vv_p_tm = conv2(vv_p_tm, kernel, 'same');    
    uu_p_tm = uu_p_tm(2:end-1, 2:end-1);
    vv_p_tm = vv_p_tm(2:end-1, 2:end-1);
    
    % get new sample grid, if resolution changes for next pass
    if (samplen(pp) ~= samplen(pp+1)) || (sampspc(pp) ~= sampspc(pp+1)) 
        [rvec, cvec] = yalebox_piv_sample_grid(samplen(pp+1), sampspc(pp+1), size(ini_ti));
        [cc_p_tm, rr_p_tm] = meshgrid(cvec, rvec);     
    end
    
    % interpolate/extrapolate to (new) predictor grid roi, restores any values lost
    % ...in smoothing , and updates the grid resolution if it has changed.
    from = ~isnan(uu_p_tm);
    to = ~from;    
    uu_p_tm(to) = spline2d(cc_p_tm(to), rr_p_tm(to), cc_p_tm(from), rr_p_tm(from), uu_p_tm(from), tension);
    vv_p_tm(to) = spline2d(cc_p_tm(to), rr_p_tm(to), cc_p_tm(from), rr_p_tm(from), vv_p_tm(from), tension);
    
    % interpolate/extrapolate displacement predictor to image resolution at
    % initial and final time
    if pp < np
        
        % propagate points to initial and final time, then re-grid to image resolution
        cc_p_ti = cc_p_tm-0.5*uu_p_tm;
        rr_p_ti = rr_p_tm-0.5*vv_p_tm;
        uu_i_ti(:) = spline2d(cc_i(:), rr_i(:), cc_p_ti(:), rr_p_ti(:), uu_p_tm(:), tension);
        vv_i_ti(:) = spline2d(cc_i(:), rr_i(:), cc_p_ti(:), rr_p_ti(:), vv_p_tm(:), tension);
        
        cc_p_tf = cc_p_tm+0.5*uu_p_tm;
        rr_p_tf = rr_p_tm+0.5*vv_p_tm;
        uu_i_tf(:) = spline2d(cc_i(:), rr_i(:), cc_p_tf(:), rr_p_tf(:), uu_p_tm(:), tension);
        vv_i_tf(:) = spline2d(cc_i(:), rr_i(:), cc_p_tf(:), rr_p_tf(:), vv_p_tm(:), tension);
    
        % deform images to midpoint time
        ini_tm = imwarp(ini_ti, -0.5*cat(3, uu_i_ti, vv_i_ti), 'cubic', 'FillValues', 0);
        fin_tm = imwarp(fin_tf,  0.5*cat(3, uu_i_tf, vv_i_tf), 'cubic', 'FillValues', 0);
        
        % deform roi masks to midpoint time, re-apply to clean up warping edge artefacts
        tmp = imwarp(double(ini_roi_ti), -0.5*cat(3, uu_i_ti, vv_i_ti), 'cubic', 'FillValues', 0);
        ini_roi_tm = abs(tmp-1) < roi_epsilon;
        ini_tm(~ini_roi_tm) = 0;
        
        tmp = imwarp(double(fin_roi_tf), 0.5*cat(3, uu_i_tf, vv_i_tf), 'cubic', 'FillValues', 0);
        fin_roi_tm = abs(tmp-1) < roi_epsilon;
        fin_tm(~fin_roi_tm) = 0;
        
    end
    
end   
% end multipass loop

uu = uu_p_tm;
vv = vv_p_tm;

% convert displacements to world coordinates 
uu = uu.*(xx(2)-xx(1)); 
vv = vv.*(yy(2)-yy(1));

% interpolate world coordinates for displacement vectors
xx = interp1(1:size(ini_ti,2), xx, cc_p_tm(1,:), 'linear', 'extrap');
yy = interp1(1:size(ini_ti,1), yy, rr_p_tm(:,1), 'linear', 'extrap');

end

%% subroutines

function [slen_ex, ilen_ex, sspc_ex] = expand_grid_def(slen, ilen, sspc, np)
%
% Expand the grid definition vectors to include the correct number of passes for
% each grid. Input arguments are defined above, but use shortened names here:
% samplen -> slen, intrlen -> ilen, sampspc -> sspc, npass -> np. 
%
% Note: outputs are intentionally not preallocated - these vectors are small and
% the performace cost is negligible. 
% %

slen_ex = [];
ilen_ex = [];
sspc_ex = [];

for ii = 1:length(np)
   slen_ex = [slen_ex, repmat(slen(ii), 1, np(ii))]; %#ok!
   ilen_ex = [ilen_ex, repmat(ilen(ii), 1, np(ii))]; %#ok!
   sspc_ex = [sspc_ex, repmat(sspc(ii), 1, np(ii))]; %#ok!
end

% repeat the last element to simplify interpolation code for the final pass
slen_ex(end+1) = slen_ex(end);
ilen_ex(end+1) = ilen_ex(end);
sspc_ex(end+1) = sspc_ex(end);

end

function [] = check_input(ini, fin, ini_roi, fin_roi, xx, yy, samplen, ...
    sampspc, intrlen, npass, valid_max, valid_eps, verbose)
%
% Check for sane input argument properties, exit with error if they do not
% match expectations.
% %

[nr, nc] = size(ini); % image size
ng = numel(samplen); % number of grid refinement steps

validateattributes(ini, {'double'}, {'2d', 'real', 'nonnan', '>=', 0, '<=' 1});
validateattributes(fin, {'double'}, {'2d', 'real', 'nonnan', '>=', 0, '<=' 1, ...
    'size', [nr, nc]});
validateattributes(ini_roi, {'logical'}, {'2d', 'size', [nr, nc]});
validateattributes(fin_roi, {'logical'}, {'2d', 'size', [nr, nc]});
validateattributes(xx, {'double'}, {'vector', 'real', 'nonnan', 'numel', nc});
validateattributes(yy, {'double'}, {'vector', 'real', 'nonnan', 'numel', nr});
validateattributes(samplen, {'numeric'}, {'vector', 'integer', 'positive', ...
    'nonnan'});
validateattributes(sampspc, {'numeric'}, {'vector', 'numel', ng, 'integer', ...
    'positive', 'nonnan'});
validateattributes(intrlen, {'numeric'}, {'vector', 'numel', ng, 'integer', ...
    'positive', 'nonnan'});
validateattributes(npass, {'numeric'}, {'vector', 'numel', ng, 'integer', ...
    'positive'});
validateattributes(valid_max, {'double'}, {'scalar', 'positive'});
validateattributes(valid_eps, {'double'}, {'scalar', 'positive'});
validateattributes(verbose, {'numeric', 'logical'}, {'scalar', 'binary'});

end

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
plot(ccnt, rcnt, 'Color', 'k', 'Marker', '.')
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
plot(ccnt, rcnt, 'Color', 'k', 'Marker', '.')
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
plot(cpk, rpk, 'Color', 'k', 'Marker', '.')
title('cross-correlation');
hold off
axis equal
axis tight

end

function [] = show_valid(drop, uu, vv)
% display results from validation step

[ii_drop, jj_drop] = find(drop);
drop_nan = ones(size(drop));
drop_nan(drop) = NaN;

subplot(2,2,1)
title('uu')
imagesc(uu)
hold on
plot(jj_drop, ii_drop, '.w')
hold off

subplot(2,2,2)
title('vv')
imagesc(vv);
hold on
plot(jj_drop, ii_drop, '.w')
hold off

subplot(2,2,3)
title('uu drop')
imagesc(uu.*drop_nan)
hold on
plot(jj_drop, ii_drop, '.w')
hold off

subplot(2,2,4)
title('vv drop')
imagesc(vv.*drop_nan);
hold on
plot(jj_drop, ii_drop, '.w')
hold off

end
