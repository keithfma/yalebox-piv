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

% variable naming scheme key
%
% ti = time of the initial image
% tm = midpoint between initial and final images
% tf = time of the final image
%
% grd = regular sample grid
% pts = irregularly spaced points
% img = regular grid at image resolution

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
[r_vec, c_vec] = yalebox_piv_sample_grid(samplen(1), sampspc(1), size(ini_ti));
[c_grd, r_grd] = meshgrid(c_vec, r_vec);
[c_img, r_img] = meshgrid(1:size(ini_ti, 2), 1:size(ini_ti ,1));

% init displacement matrices
sz = size(c_grd);
u_grd_tm = zeros(sz); 
v_grd_tm = zeros(sz); 
du_pts_tm = zeros(sz); 
dv_pts_tm = zeros(sz); 

sz = size(ini_ti);
u_img_ti = zeros(sz);
v_img_ti = zeros(sz);
u_img_tf = zeros(sz);
v_img_tf = zeros(sz);

% init deformed images
ini_tm = ini_ti;
fin_tm = fin_tf;

% multipass loop
np = length(samplen); 
for pp = 1:np
    
    % reset per-pass variables
    sz = size(c_grd);
    c_pts = zeros(sz);
    r_pts = zeros(sz);    
    du_pts_tm = nan(sz);
    dv_pts_tm = nan(sz);
    du_grd_tm = nan(sz);
    dv_grd_tm = nan(sz);    
    roi = true(sz);
    min_overlap = min_frac_overlap*samplen(pp)*samplen(pp);    
    
    % get corrector displacements on the predictor grid
    for kk = 1:numel(du_pts_tm)
            
        % get sample window and its centroid location
        [samp, samp_pos, frac_data, r_samp_cntr, c_samp_cntr] = ...
            yalebox_piv_window(ini_tm, r_grd(kk), c_grd(kk), samplen(pp));
        
        % get interrogation window
        [intr, intr_pos] = ...
            yalebox_piv_window(fin_tm, r_grd(kk), c_grd(kk), intrlen(pp));
        
        % skip if sample window is too empty
        if frac_data < min_frac_data
            roi(kk) = false;
            du_pts_tm(kk) = NaN;
            dv_pts_tm(kk) = NaN;
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
        du_pts_tm(kk) = cpeak-samplen(pp)-(samp_pos(1)-intr_pos(1));
        dv_pts_tm(kk) = rpeak-samplen(pp)-(samp_pos(2)-intr_pos(2));
        
        % compute location of this observation at midpoint time
        c_pts(kk) = c_samp_cntr + 0.5*du_pts_tm(kk);
        r_pts(kk) = r_samp_cntr + 0.5*dv_pts_tm(kk);
                
    end 
    
    % validate dislacement updates 
    [du_pts_tm, dv_pts_tm] = piv_validate_pts_nmed(c_pts, r_pts, du_pts_tm, dv_pts_tm, 8, valid_max, valid_eps, true);    
    
    % interpolate displacement update to sample grid, points outside roi remain NaN
    from = ~isnan(du_pts_tm); 
    to = roi;
    % du_grd_tm(to) = spline2d(c_grd(to), r_grd(to), c_pts(from), r_pts(from), du_pts_tm(from), tension);
    % dv_grd_tm(to) = spline2d(c_grd(to), r_grd(to), c_pts(from), r_pts(from), dv_pts_tm(from), tension);
    
    span_pts = 9;
    span_frac = span_pts/sum(from(:));
    
    % % dataviz version: many matrix scaling warnings, suggests a poor implementation
    % du_grd_tm(:) = NaN;
    % du_grd_tm(to) = loess2(c_pts(from), r_pts(from), du_pts_tm(from), c_grd(to), r_grd(to), span_frac, 1, true);
    % dv_grd_tm(:) = NaN;
    % dv_grd_tm(to) = loess2(c_pts(from), r_pts(from), dv_pts_tm(from), c_grd(to), r_grd(to), span_frac, 1, true);
    
    % matlab version
    model = fit([c_pts(from), r_pts(from)], du_pts_tm(from), 'lowess', 'Span', span_frac, 'Robust', 'bisquare');    
    du_grd_tm(:) = NaN;
    du_grd_tm(to) = model(c_grd(to), r_grd(to));
    model = fit([c_pts(from), r_pts(from)], dv_pts_tm(from), 'lowess', 'Span', span_frac, 'Robust', 'bisquare');    
    dv_grd_tm(:) = NaN;
    dv_grd_tm(to) = model(c_grd(to), r_grd(to));

    keyboard
    
    % update displacement, points outside roi become NaN
    u_grd_tm = u_grd_tm + du_grd_tm;
    v_grd_tm = v_grd_tm + dv_grd_tm;
    
    % NOTE: could make better use of the edge data by accounting for
    % displacement by the smoothing kernel.
    
%     % smooth predictors, 3x3 kernel smoother...
%     % % NaNs at all roi boundaries propagate inward to all points affected by
%     % % the bounday
%     u_grd_tm = padarray(u_grd_tm, [1 1], NaN, 'both');    
%     v_grd_tm = padarray(v_grd_tm, [1 1], NaN, 'both');        
%     kernel = fspecial('average', 3);    
%     u_grd_tm = conv2(u_grd_tm, kernel, 'same');
%     v_grd_tm = conv2(v_grd_tm, kernel, 'same');    
%     u_grd_tm = u_grd_tm(2:end-1, 2:end-1);
%     v_grd_tm = v_grd_tm(2:end-1, 2:end-1);
    
    % interpolate/extrapolate points lost in smoothing
    from = ~isnan(u_grd_tm);
    to = roi;    
    u_grd_tm(to) = spline2d(c_grd(to), r_grd(to), c_grd(from), r_grd(from), u_grd_tm(from), tension);
    v_grd_tm(to) = spline2d(c_grd(to), r_grd(to), c_grd(from), r_grd(from), v_grd_tm(from), tension);

    % prepare for next pass
    if pp < np
        
        % extend sample grid to cover full image footprint, no need to be regular
        c_vec_ext = [1; c_vec(:); c_img(1,end)];
        r_vec_ext = [1; r_vec(:); r_img(end,1)];
        [c_ext, r_ext] = meshgrid(c_vec_ext, r_vec_ext);
        
        % populate extended sample grid using tension splines        
        from = ~isnan(u_grd_tm);
        u_ext_tm = zeros(size(c_ext));
        v_ext_tm = zeros(size(c_ext));
        u_ext_tm(:) = spline2d(c_ext(:), r_ext(:), c_grd(from), r_grd(from), u_grd_tm(from), tension);
        v_ext_tm(:) = spline2d(c_ext(:), r_ext(:), c_grd(from), r_grd(from), v_grd_tm(from), tension);

        % re-grid to image resolution at initial time using cheaper interpolant
        c_pts = c_ext - 0.5*u_ext_tm;
        r_pts = r_ext - 0.5*v_ext_tm;        
        
        interpolant = scatteredInterpolant(c_pts(:), r_pts(:), u_ext_tm(:), 'natural', 'linear');
        u_img_ti = interpolant(c_img, r_img);
        interpolant.Values = v_ext_tm(:);
        v_img_ti = interpolant(c_img, r_img);
        
        % re-grid to image resolution at initial time using cheaper interpolant
        c_pts = c_ext + 0.5*u_ext_tm;
        r_pts = r_ext + 0.5*v_ext_tm;        
        
        interpolant = scatteredInterpolant(c_pts(:), r_pts(:), u_ext_tm(:), 'natural', 'linear');
        u_img_tf = interpolant(c_img, r_img);
        interpolant.Values = v_ext_tm(:);
        v_img_tf = interpolant(c_img, r_img);
        
        % deform images to midpoint time
        ini_tm = imwarp(ini_ti, -0.5*cat(3, u_img_ti, v_img_ti), 'cubic', 'FillValues', 0);
        fin_tm = imwarp(fin_tf,  0.5*cat(3, u_img_tf, v_img_tf), 'cubic', 'FillValues', 0);
        
        % deform roi masks to midpoint time, re-apply to clean up warping edge artefacts
        tmp = imwarp(double(ini_roi_ti), -0.5*cat(3, u_img_ti, v_img_ti), 'cubic', 'FillValues', 0);
        ini_roi_tm = abs(tmp-1) < roi_epsilon;
        ini_tm(~ini_roi_tm) = 0;
        
        tmp = imwarp(double(fin_roi_tf), 0.5*cat(3, u_img_tf, v_img_tf), 'cubic', 'FillValues', 0);
        fin_roi_tm = abs(tmp-1) < roi_epsilon;
        fin_tm(~fin_roi_tm) = 0;
        
        % part 2: get new coordinate and predictor grids if resolution changes
        % in next pass
        
        if (samplen(pp) ~= samplen(pp+1)) || (sampspc(pp) ~= sampspc(pp+1))
            
            [r_vec, c_vec] = yalebox_piv_sample_grid(samplen(pp+1), sampspc(pp+1), size(ini_ti));
            [c_grd, r_grd] = meshgrid(c_vec, r_vec);
            
            fprintf('F\n');
            from = ~isnan(u_grd_tm);
            u_grd_tm = spline2d(c_grd(:), r_grd(:), c_grd(from), r_grd(from), u_grd_tm(from), tension);
            v_grd_tm = spline2d(c_grd(:), r_grd(:), c_grd(from), r_grd(from), v_grd_tm(from), tension);
            
            u_grd_tm = reshape(u_grd_tm, size(c_grd));
            v_grd_tm = reshape(v_grd_tm, size(c_grd));           
        end        
    end
    
end   
% end multipass loop

% rename displacements
uu = u_grd_tm;
vv = v_grd_tm;

% convert displacements to world coordinates 
uu = uu.*(xx(2)-xx(1)); 
vv = vv.*(yy(2)-yy(1));

% interpolate world coordinates for displacement vectors
xx = interp1(1:size(ini_ti,2), xx, c_grd(1,:), 'linear', 'extrap');
yy = interp1(1:size(ini_ti,1), yy, r_grd(:,1), 'linear', 'extrap');

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
