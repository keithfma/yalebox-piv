function out = post_strain(x, y, uu, vv, roi, pad_method)
% function out = post_strain(x, y, uu, vv, roi, pad_method)
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
%   out: Struct containing results, attribures follow.
%   out.F11: Deformation-gradient tensor, element (1,1), dx/dX, [1]
%   out.F12: Deformation-gradient tensor, element (1,2), dx/dY, [1]
%   out.F21: Deformation-gradient tensor, element (2,1), dy/dX, [1]
%   out.F22: Deformation-gradient tensor, element (2,2), dy/dY, [1]
%   out.S1: Maximum principle stretch, [1]
%   out.S2: Minimum principle stretch, [1]
%   out.S1x: x-component of maximum principle stretch direction
%   out.S1y: y-component of maximum principle stretch direction
%   out.S2x: x-component of minimum principle stretch direction
%   out.S2y: y-component of minimum principle stretch direction
%   out.spin: Rotation angle, clockwise, [radians]
%   out.D1: Maximum equivalent stretching rate, [1/step]
%   out.D2: Minimum equivalent stretching rate, [1/step]
%   out.Dt: Scalar total strain rate, as in (refs 1-3), [1/step]
%   out.Dd: Scalar deviatoric strain rate, as in  (refs 1-3), [1/step]
%   out.Dv: Volume strain rate, as in (refs 1-3), [1/step]
%   out.Wk: Kinematic vorticity number relative to Dt, as in refs (1-3), [1]
%   out.Wk_star: " " relative to Dd, as in refs (1-3), [1]
%   out.Ak: Kinematic dilatancy number relative to Dt, as in refs (1-3), [1]
%   out.Ak_star: " " relative to Dd, as in refs (1-3), [1]
%
% References:
%   (1) Brandon, M. T. (2015). Instantaneous and Finite Descriptions of a
%       Velocity Field. Unpublished Technical Report.
%   (2) Brandon, M. T. (1995). Analysis of geologic strain data in
%       strain-magnitude space. Journal of Structural Geology, 17(10),
%       1375–1385. http://doi.org/10.1016/0191-8141(95)00032-9
%   (3) Ring, U., & Brandon, M. T. (1999). Ductile deformation and mass
%       loss in the Franciscan Subduction Complex: implications for
%       exhumation processes in accretionary wedges. Exhumation Processes:
%       Normal Faulting, Ductile Flow and Erosion, (154) 55-86.
%       http://doi.org/10.1144/gsl.sp.1999.154.01.03
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

% interpolate displacement vectors at initial time 
%   this means the x,y grid represents initial position (i.e. X)
[xx, yy] = meshgrid(x, y);
xi_pts = xx(roi)-uu(roi)/2;
yi_pts = yy(roi)-vv(roi)/2;
tension = 0.9;
uu(roi) = spline2d(xx(roi), yy(roi), xi_pts, yi_pts, uu(roi), tension);
vv(roi) = spline2d(xx(roi), yy(roi), xi_pts, yi_pts, vv(roi), tension);

% compute deformation gradient tensor: F = dx/dX = d/dX( X+u ) = I + du/dX
[F11, F12] = spatial_gradient(x, y, uu, pad_method);
F11 = F11 + 1;
[F21, F22] = spatial_gradient(x, y, vv, pad_method);
F22 = F22+1;

% preallocate derived strain parameters
S1 = nan(ny, nx);
S2 = nan(ny, nx);
S1x = nan(ny, nx);
S1y = nan(ny, nx);
S2x = nan(ny, nx);
S2y = nan(ny, nx);
spin = nan(ny, nx);

% calculate derived strain parameters
for kk = find(roi)'
    F = [F11(kk), F12(kk); F21(kk), F22(kk)];
    % polar decomposition, F = VR, B = V^2 = FF', and R = (V^-1)F,
    % %     where B is the Left Cauchy-Green tensor.
    B = F*F';
    % % eigen solution sorted with eigenvalues in descending order
    [T, lambda] = eig(B);
    [lambda, order] = sort(diag(lambda),'descend');
    S1(kk) = sqrt(lambda(1));
    S2(kk) = sqrt(lambda(2));
    T = T(:, order);
    % store unit vectors for maximum and minimum extension rate direction
    S1x(kk) = T(1,1);
    S1y(kk) = T(2,1);
    S2x(kk) = T(1,2);
    S2y(kk) = T(2,2);
    % rotation matrix and spin
    V_inverse = T*diag(lambda.^(-1/2))*T';
    R = V_inverse*F;
    spin(kk) = atan2(R(2,1), R(1,1));
end
% equivalent steady rate-of-deformation / stretching tensor components
%   NOTE: D = logm(S), which is trivial for S in principle form
D1 = log(S1);
D2 = log(S2);

% calculate invariants, assuming plane strain (D3=0)
Dt = sqrt(D1.^2 + D2.^2);
Dd = sqrt(((D1-D2).^2 + D1.^2 + D2.^2)/3);
Dv = sqrt(1/3)*(D1 + D2); 
 
% calculate kinematic numbers
Wk      = 2*spin./(sqrt(2)*Dt);
Wk_star = 2*spin./(sqrt(2)*Dd);
Ak      = Dv./(sqrt(2)*Dt);
Ak_star = Dv./(sqrt(2)*Dd);
% NOTE: div-by-zero creates inf values that cannot be written to netCDF  
Wk(isinf(Wk)) = NaN;
Wk_star(isinf(Wk_star)) = NaN;
Ak(isinf(Ak)) = NaN;
Ak_star(isinf(Ak_star)) = NaN;

% prepare outputs
out = struct();
out.F11 = F11;
out.F12 = F12;
out.F21 = F21;
out.F22 = F22;
out.S1 = S1;
out.S2 = S2;
out.S1x = S1x;
out.S1y = S1y;
out.S2x = S2x;
out.S2y = S2y;
out.spin = spin;
out.D1 = D1;
out.D2 = D2;
out.Dt = Dt;
out.Dd = Dd;
out.Dv = Dv;
out.Wk = Wk;
out.Ak = Ak;
out.Wk_star = Wk_star;
out.Ak_star = Ak_star;

function [dzdx, dzdy] = spatial_gradient(x, y, zz, pad_method)
%
% Return numerical gradient using optimal 7-tap method from reference (1).
%
% Arguments:
%   x, y: Coordinate vectors for x- and y-directions
%   zz: Data matrix for which gradient is to be computed, assumes missing
%       data (outside the ROI) is NaN
%   pad_method: Select padding method, valid options are {'nearest'}
%   dzdx, dzdy: x- and y-direction components of the gradient of zz
%
% References:
%   (1) Farid, H., & Simoncelli, E. P. (2004). Differentiation of discrete
%       multidimensional signals. IEEE Transactions on Image Processing,
%       13(4), 496–508. http://doi.org/10.1109/TIP.2004.823819
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