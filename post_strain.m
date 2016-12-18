function [L, F, S] = post_strain(x, y, uu, vv, roi, pad_method)
%
% Compute strain deformation parameters for input velocity fields.
%
% NOTE: Assumes x and y are regularly spaced (i.e. dx and dy are constant)
%
% Arguments:
%   x, y: Coordinate vectors for x- and y-directions
%   uu, vv: Displacement field matrices for x- and y-direction, [m/step]
%   roi: Region-of-interest mask matrix (false for no data)
%   pad_method: (temporary) Select padding method, valid options are 
%       {'nearest'}
%   L: Displacement-gradient tensor field, 3D matrix, with the 3rd
%       dimension containing tensor components ordered column-wise, e.g.
%       (L11, L21, L12, L22), [1]
%   F: Deformation-gradient tensor, 3D matrix, stored same as L, [1]
%   S: Principle stretches, 3D matrix, S(:,:,1) contains the maximum
%       principle stretch, [1]
%   D: Principle natural strains (log(S)), 3D matrix, ordered as in S, [1]
% %

% check inputs
narginchk(6,7);
validateattributes(x, {'numeric'}, {'vector','real'}, mfilename, 'x');
validateattributes(y, {'numeric'}, {'vector','real'}, mfilename, 'y');
nx = length(x);
ny = length(y);
validateattributes(uu, {'numeric'}, {'size', [ny, nx]}, mfilename, 'uu');
validateattributes(vv, {'numeric'}, {'size', [ny, nx]}, mfilename, 'vv');
validateattributes(roi,  {'logical'}, {'size', [ny, nx]}, mfilename, 'roi');
validateattributes(pad_method, {'char'}, {'vector'}, mfilename, 'pad_method');

% <DEBUG> compute F directly
% consider coordinate grid (x, y) to be initial position, get final
%   position for points on this grid
[XX, YY] = meshgrid(x, y);
xi_pts = XX(roi)-uu(roi)/2;
yi_pts = YY(roi)-vv(roi)/2;
xf_pts = XX(roi)+uu(roi)/2;
yf_pts = YY(roi)+vv(roi)/2;
tension = 0.9;
xx = nan(ny, nx);
yy = nan(ny, nx);
xx(roi) = spline2d(XX(roi), YY(roi), xi_pts, yi_pts, xf_pts, tension);
yy(roi) = spline2d(XX(roi), YY(roi), xi_pts, yi_pts, yf_pts, tension);

[f11, f12] = spatial_gradient(x, y, xx, pad_method);
[f21, f22] = spatial_gradient(x, y, yy, pad_method);
% </DEBUG>

% calculate displacement-gradient tensor 
[L11, L12] = spatial_gradient(x, y, uu, pad_method);
[L21, L22] = spatial_gradient(x, y, vv, pad_method);

% calculate deformation-gradient tensor 
% NOTE: F = L + I, where L is the displacement-gradient tensor 
F11 = L11 + 1;
F12 = L12;
F21 = L21;
F22 = L22 + 1;

% NOTE: compared two versions of F, they are basically the same except
% slightly shifted, there is some funny business at the edge...


keyboard


% preallocate derived strain parameters
S1 = nan(ny, nx);
S2 = nan(ny, nx);
S1x = nan(ny, nx);
S1y = nan(ny, nx);
S2x = nan(ny, nx);
S2y = nan(ny, nx);

% calculate derived strain parameters
for kk = find(roi)'
    F = [F11(kk), F12(kk); F21(kk), F22(kk)];
    % % polar decomposition, F = VR, B = V^2 = FF', and R = (V^-1)F,
    % %     where B is the Left Cauchy-Green tensor.
    B = F*F';
    % % eigen solution sorted with eigenvalues in descending order
    [T, lambda] = eig(B);
    [lambda, order] = sort(diag(lambda),'descend');
    S1(kk) = sqrt(lambda(1));
    S2(kk) = sqrt(lambda(2));
    T = T(:, order);
    % % store unit vectors for maximum and minimum extension rate direction
    S1x(kk) = T(1,1);
    S1y(kk) = T(2,1);
    S2x(kk) = T(1,2);
    S2y(kk) = T(2,2);
    % % rotation matrix and spin
end
D1 = log(S1);
D2 = log(S2);

% %... Calculate the spin tensor W, and the stretch tensor D in its principal form
% % D1 and D2 are the principal extension and shortening rates, respectively
% % theta contains the principal shortening rate directions in radians
% D1 = nan(mPad,nPad);
% D2 = nan(mPad,nPad);
% D2x = nan(mPad,nPad);
% D2y = nan(mPad,nPad);
% spin = nan(mPad,nPad);
% for i = 1:mPad
%     for j = 1:nPad
%         %... Skip if F contains one or more nans
%         if ~isnan(F11(i,j)*F12(i,j)*F21(i,j)*F22(i,j))
%             F = [F11(i,j), F12(i,j); F21(i,j),  F22(i,j)];
%             %... Polar decomposition, F = VR, B = V^2 = FF', and R = (V^-1)F,
%             % where B is the Left Cauchy-Green tensor.
%             B = F*F';
%             %... Calculate eigen solution and sort with eigenvalues in descending order
%             [T,lambda] = eig(B);
%             [lambda, order] = sort(diag(lambda),'descend');
%             S = sqrt(lambda);
%             D1(i,j) = log(S(1));
%             D2(i,j) = log(S(2));
%             T = T(:,order);
%             %... Calculate unit vector for maximum extension rate direction
%             theta = acos(T(1,1));
%             D2x(i,j) = cos(theta);
%             D2y(i,j) = sin(theta);
%             %... Calculate spin
%             invV = T*diag(lambda.^(-1/2))*T';
%             R = invV*F;
%             spin(i,j) = atan2(R(2,1),R(1,1));
%         end
%     end
% end
% clear F11 F12 F21 F22
% 
% %... Calculate invariants, assuming plane strain (D3=0)
% % Dt = sqrt(D1.^2 + D2.^2);
% Dd = sqrt(((D1-D2).^2 + D1.^2 + D2.^2)/3);
% Dv = sqrt(1/3)*(D1 + D2);
% 
% %... Smooth displacement field using Gaussian kernel with sigma = 1.2
% % pixels, as recommended on p. 512 in Adrian and Westerweel, 2011,
% % Particle Image Velocimetry. The filter in imgaussian defaults
% % to a size equal to +/-3sigma.
% % spin = imgaussian(spin,1.2,3);
% % Dv = imgaussian(Dv,1.2,3);
% % Dd = imgaussian(Dd,1.2,3);
% % % Dt = imgaussian(Dt,1.2,3);
% % D2x = imgaussian(D2x,1.2,3);
% % D2y = imgaussian(D2y,1.2,3);
% 
% %... Remove padding
% uX = uX((1:m)+padSize,:);
% uY = uY((1:m)+padSize,:);
% spin = spin((1:m)+padSize,:);
% Dv = Dv((1:m)+padSize,:);
% Dd = Dd((1:m)+padSize,:);
% % Dt = Dt((1:m)+padSize,:);
% D2x = D2x((1:m)+padSize,:);
% D2y = D2y((1:m)+padSize,:);
% uX(~iUXY) = nan;
% uY(~iUXY) = nan;
% spin(~iUXY) = nan;
% Dv(~iUXY) = nan;
% Dd(~iUXY) = nan;
% % Dt(~iUXY) = nan;
% D2x(~iUXY) = nan;
% D2y(~iUXY) = nan;
% 
% %... Calculate displacement magnitude
% displacement = sqrt(uX.^2 + uY.^2);
% 
% %... Calculate kinematic numbers
% % Wk = 2*spin./(sqrt(2)*Dt);
% % Ak = Dv./(sqrt(2)*Dt);
% WkStar = 2*spin./(sqrt(2)*Dd);
% AkStar = Dv./(sqrt(2)*Dd);



% prepare outputs
L = cat(3, L11, L21, L12, L22); % tensor store comlum-wise in 3rd dim
F = cat(3, F11, F21, F12, F22);
S = cat(3, S1, S2);
D = cat(3, D1, D2);

function [dzdx, dzdy] = spatial_gradient(x, y, zz, pad_method)
%
% Return numerical gradient using optimal 7-tap method from reference [1].
%
% Arguments:
%   x, y: Coordinate vectors for x- and y-directions
%   zz: Data matrix for which gradient is to be computed, assumes missing
%       data (outside the ROI) is NaN
%   pad_method: Select padding method, valid options are {'nearest'}
%   dzdx, dzdy: x- and y-direction components of the gradient of zz
%
% References:
% [1] 
% 
% Dependencies:
%   derivativesByFilter 
% %

% define parameters
pad_width = 3;

% get some constants
dx = x(2)-x(1);
dy = y(2)-y(1);
roi = ~isnan(zz);

% pad coordinate vectors
pad_coord = @(z, dz) [z(1)+dz*(-pad_width:-1)'; z(:); z(end)+dz*(1:pad_width)'];
x_p = pad_coord(x, dx);
y_p = pad_coord(y, dy);
[xx_p, yy_p] = meshgrid(x_p, y_p);

% pad data matrix
% TODO: add Brandon padding
% TODO: add inpaint_nans padding
zz_p = padarray(zz, [pad_width, pad_width], NaN, 'both');
roi_p = padarray(roi, [pad_width, pad_width], false, 'both');

if strcmp(pad_method, 'nearest')
    si = scatteredInterpolant(xx_p(roi_p), yy_p(roi_p), zz_p(roi_p), ...
        'nearest', 'nearest');
    zz_p(~roi_p) = si(xx_p(~roi_p), yy_p(~roi_p));
else
    error('invalid padding method selected');
end

% compute gradient
[dzdx_p, dzdy_p] = derivativesByFilters(zz_p, 'x', 'y', dx, dy, 'seven');

% remove pad
unpad = @(m) m(pad_width+1:end-pad_width, pad_width+1:end-pad_width);
dzdx = unpad(dzdx_p);
dzdy = unpad(dzdy_p);

% re-apply ROI mask
dzdx(~roi) = NaN;
dzdy(~roi) = NaN;