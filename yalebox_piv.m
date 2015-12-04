function [xx, yy, uu, vv] = ...
    yalebox_piv(ini, fin, xx, yy, samplen, sampspc, intrlen, ...
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

% temporary
warning off images:removing:function

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
rr_full = 1:size(ini, 1);
cc_full = 1:size(ini, 2);

% init sample grid 
[rr, cc] = yalebox_piv_sample_grid(samplen(1), sampspc(1), size(ini));
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
        print_sep(sprintf('grid refinement step %i of %i', gg, ngrid));
    end
    
    % copy old sample grid
    rr_old = rr;
    cc_old = cc;
    
    % get new sample grid    
    [rr, cc] = yalebox_piv_sample_grid(samplen(gg), sampspc(gg), size(ini));
    nr = length(rr);
    nc = length(cc);
    
    % interpolate/extrapolate displacements from old to new sample grid    
    % CHECK HERE
    [uu, vv] = yalebox_piv_interp2d(rr_old, cc_old, uu, vv, rr, cc, 'spline');
    
    % loop over image deformation passes
    for pp = 1:npass(gg)
        
        % report image deformation step
        if verbose
            print_sep(sprintf('\timage deformation pass %i of %i', pp, npass(gg)));
        end
        
        % interpolate/extrapolate displacement vectors to full image resolution
        % CHECK HERE
        [uu_full, vv_full] = yalebox_piv_interp2d(rr, cc, uu, vv, rr_full, cc_full, 'spline');
        
        % deform images (does nothing if uu0 and vv0 are 0)
        % CHECK HERE: deformation does not yield clean edges
        defm_ini = imwarp(ini, -cat(3, uu_full, vv_full)/2, 'cubic', 'FillValues', 0);
        defm_fin = imwarp(fin,  cat(3, uu_full, vv_full)/2, 'cubic', 'FillValues', 0);
               
        % all grid points start in the ROI
        roi = true(nr, nc);
        
        % set subpixel correlation value matrix to zero
        cval = zeros(nr, nc);
        
        % reset data centroid grids
        rr_cntr = zeros(nr, nc);
        cc_cntr = zeros(nr, nc);
        
        % loop over sample grid
        for jj = 1:nc
            for ii = 1:nr
                
                % get sample and (offset) interrogation windows
                [samp, samp_pos, frac_data, rr_cntr(ii,jj), cc_cntr(ii,jj)] = ...
                    yalebox_piv_window(defm_ini, rr(ii), cc(jj), samplen(gg));
                [intr, intr_pos] = ...
                    yalebox_piv_window(defm_fin, rr(ii), cc(jj), intrlen(gg));
                   
                % skip and remove from ROI if sample window is too empty
                if frac_data < 0.25
                    roi(ii, jj) = false;
                    continue
                end                    
                
                % compute normalized cross-correlation
                xcr = normxcorr2(samp, intr);
                                
                % find correlation plane max, subpixel precision
                [rpeak, cpeak, val, stat] = yalebox_piv_peak_gauss2d(xcr);
                % [rpeak, cpeak, stat] = peak_optim_fourier(xcr);
                if stat == false
                    uu(ii, jj) = NaN;
                    vv(ii, jj) = NaN;
                    continue
                end
                
                % find displacement from position of the correlation max
                %   - account for padding in cross-correlation (-samplen(gg))
                %   - account for relative position of interogation and sample
                %     windows (e,g, for columns: -(samp_pos(1)-intr_pos(1))
                delta_uu = cpeak-samplen(gg)-(samp_pos(1)-intr_pos(1));
                delta_vv = rpeak-samplen(gg)-(samp_pos(2)-intr_pos(2));
                
                uu(ii, jj) = uu(ii, jj)+delta_uu;
                vv(ii, jj) = vv(ii, jj)+delta_vv;
                cval(ii, jj) = val;
                
                % % debug {
                % fprintf('CVAL = %f\n', cval(ii, jj));
                % figure(1)
                % show_win(defm_ini, defm_fin, rr(ii), cc(jj), samp, samp_pos, intr, intr_pos);
                % figure(2)
                % show_xcor(xcr, rpeak, cpeak);
                % pause
                % % } debug
                
            end % ii
        end % jj
        
        % find and drop invalid displacement vectors
        drop = yalebox_piv_valid_nmed(uu, vv, valid_max, valid_eps);  
        uu(drop) = NaN;
        vv(drop) = NaN;
        
        % interpolate/extrapolate from partial centroid grid to full regular grid        
        
%         % option 1: thin plate splines
%         p_smooth = 0.25; % smoothing parameter
%         [cc_grid, rr_grid] = meshgrid(cc, rr);
%         
%         warning('off', 'SPLINES:TPAPS:NaNs');        
%         st = tpaps([cc_cntr(:)'; rr_cntr(:)'], uu(:)', p_smooth);        
%         uu = reshape( fnval(st, [cc_grid(:)'; rr_grid(:)']), nr, nc);
%          
%         st = tpaps([cc_cntr(:)'; rr_cntr(:)'], vv(:)', p_smooth);        
%         vv = reshape( fnval(st, [cc_grid(:)'; rr_grid(:)']), nr, nc);        
%         warning('on', 'SPLINES:TPAPS:NaNs');
        
        % option 2: scattered interpolant
        
%         try
%             [cc_grid, rr_grid] = meshgrid(cc, rr);
%             
%             interpolant = scatteredInterpolant(cc_cntr(:), rr_cntr(:), uu(:), ...
%                 'nearest', 'nearest');
%             uu = interpolant(cc_grid, rr_grid);
%             
%             interpolant.Values = vv(:);
%             vv = interpolant(cc_grid, rr_grid);
%         catch err
%             fprintf(getReport(err));
%         end
            
        
        % debug: plot centroids and regular grid {
        imagesc(ini); colormap('gray');
        hold on
        plot(cc_cntr(:), rr_cntr(:), 'xb')
        [cc_grid, rr_grid] = meshgrid(cc, rr);
        plot(cc_grid(:), rr_grid(:), 'or')
        % } debug
        
        keyboard
        
        % % debug {
        % show_valid(drop, uu, vv);
        % pause
        % % } debug
       
    end % pp
    
end % gg

% delete points outside the ROI
uu(~roi) = NaN;
vv(~roi) = NaN;

% convert displacements to world coordinates (assumes constant grid spacing)
uu = uu.*(xx(2)-xx(1));
vv = vv.*(yy(2)-yy(1));

% interpolate world coordinates for displacement vectors
xx = interp1(1:size(ini,2), xx, cc, 'linear', 'extrap');
yy = interp1(1:size(ini,1), yy, rr, 'linear', 'extrap');

% % debug {
% keyboard
% % } debug

end

%% subroutines

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
