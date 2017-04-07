function [r1, c1, u1, v1, x1, y1] = ...
    piv_interp_sample_grid(r0, c0, u0, v0, roi, len, spc, xw, yw, tension, verbose)
% function [r1, c1, u1, v1, x1, y1] = ...
%     piv_interp_sample_grid(r0, c0, u0, v0, roi, len, spc, xw, yw, tension, verbose)
%
% Get new sample grid and interpolate data to this new grid. Interpolation fills
% the entire new grid. For additinal piv passes, the new roi is imposed later by
% adding NaNs to the current estimate.
%
% Arguments:
%
%   r0, c0 = 2D matrix, original sample grid in pixel coordinates
%
%   u0, v0 = 2D matrix, displacement estimates on original sample grid
%
%   roi = 2D matrix, logical, flag indicating if point is within the data ROI
%       (1) or not (0)
%
%   len = Scalar, integer, side length for sample windows 
%
%   spc = Scalar, integer, spacing between sample windows
%   
%   xw, yw = Vectors, world coordinates at full image resolution
%
%   tension = Scalar, tension parameter for spline interpolation.
%
%   verbose = Logical flag, display verbose messages (1) or don't
% %

% TODO: This function is obsolete, I handle call piv_interp* directly from piv()
    
[r1, c1, x1, y1] = piv_sample_grid(len, spc, xw, yw);

if verbose
    fprintf('%s: new grid = [%d, %d], old grid = [%d, %d]\n', ...
        mfilename, size(r1,1), size(r1,2), size(r0,1), size(r0,2));
end

% <DEBUG>: selected biharmonic
u1 = griddata(c0(roi), r0(roi), u0(roi), c1, r1, 'v4');
v1 = griddata(c0(roi), r0(roi), v0(roi), c1, r1, 'v4');
% </DEBUG>

% % <DEBUG: spline-in-tension>
% u1_t = spline2d(c1(:), r1(:), c0(roi), r0(roi), u0(roi), tension);
% v1_t = spline2d(c1(:), r1(:), c0(roi), r0(roi), v0(roi), tension);
% 
% u1_t = reshape(u1, size(r1));
% v1_t = reshape(v1, size(r1));
% % </DEBUG>
% 
% % <DEBUG: biharmonic>
% u1_b = griddata(c0(roi), r0(roi), u0(roi), c1, r1, 'v4');
% v1_b = griddata(c0(roi), r0(roi), v0(roi), c1, r1, 'v4');
% % </DEBUG>
% 
% % <DEBUG: cubic> 
% u1_c = griddata(c0(roi), r0(roi), u0(roi), c1, r1, 'cubic');
% v1_c = griddata(c0(roi), r0(roi), v0(roi), c1, r1, 'cubic');
% % </DEBUG>
% 
% % <DEBUG: natural> 
% u1_n = griddata(c0(roi), r0(roi), u0(roi), c1, r1, 'natural');
% v1_n = griddata(c0(roi), r0(roi), v0(roi), c1, r1, 'natural');
% % </DEBUG>
% 
% % <DEBUG: linear> 
% u1_l = griddata(c0(roi), r0(roi), u0(roi), c1, r1, 'linear');
% v1_l = griddata(c0(roi), r0(roi), v0(roi), c1, r1, 'linear');
% % </DEBUG>
% 
% % <DEBUG: comparo>
% subplot(3,2,1); 
% imagesc(u1_t)
% title('original')
% 
% subplot(3,2,2); 
% imagesc(u1_b)
% title('biharmonic')
% 
% subplot(3,2,3); 
% imagesc(u1_c)
% title('cubic')
% 
% subplot(3,2,4); 
% imagesc(u1_n)
% title('natural')
% 
% subplot(3,2,5); 
% imagesc(u1_l)
% title('linear')
% %</DEBUG>