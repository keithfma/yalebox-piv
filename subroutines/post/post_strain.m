function out = post_strain(x, y, uu, vv)
% function out = post_strain(x, y, uu, vv)
%
% Compute strain parameters for input velocity fields.
%
% NOTE: Assumes x and y are regularly spaced (i.e. dx and dy are constant)
%
% Arguments:
%   x, y: Coordinate vectors for x- and y-directions
%   uu, vv: Displacement field matrices for x- and y-direction, [m/step]
%   out: Struct containing results, attributes follow.
%   out.x, out.y: See x, y above
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
%       1375?1385. http://doi.org/10.1016/0191-8141(95)00032-9
%   (3) Ring, U., & Brandon, M. T. (1999). Ductile deformation and mass
%       loss in the Franciscan Subduction Complex: implications for
%       exhumation processes in accretionary wedges. Exhumation Processes:
%       Normal Faulting, Ductile Flow and Erosion, (154) 55-86.
%       http://doi.org/10.1144/gsl.sp.1999.154.01.03
% %

% add dependencies
update_path('spline', 'deriv');

% check inputs
narginchk(4, 4);
validateattributes(x, {'numeric'}, {'vector','real'}, mfilename, 'x');
validateattributes(y, {'numeric'}, {'vector','real'}, mfilename, 'y');
nx = length(x);
ny = length(y);
validateattributes(uu, {'numeric'}, {'size', [ny, nx]}, mfilename, 'uu');
validateattributes(vv, {'numeric'}, {'size', [ny, nx]}, mfilename, 'vv');

% compute deformation gradient tensor: F = dx/dX = d/dX( X+u ) = I + du/dX
[F11, F12] = spatial_gradient(x, y, uu);
F11 = F11 + 1;
[F21, F22] = spatial_gradient(x, y, vv);
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
for kk = 1:numel(uu)
    
    % nothing more to do if strain matrix has nans
    if any(isnan([F11(kk), F12(kk), F21(kk), F22(kk)]))
        continue
    end
    
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
out.x = x;
out.y = y;
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

function [dzdx, dzdy] = spatial_gradient(x, y, zz)
%
% Return numerical gradient using optimal 7-tap method from reference (1).
%
% Arguments:
%   x, y: Coordinate vectors for x- and y-directions
%   zz: Data matrix for which gradient is to be computed, assumes missing
%       data (outside the ROI) is NaN
%   dzdx, dzdy: x- and y-direction components of the gradient of zz
%
% References:
%   (1) Farid, H., & Simoncelli, E. P. (2004). Differentiation of discrete
%       multidimensional signals. IEEE Transactions on Image Processing,
%       13(4), 496?508. http://doi.org/10.1109/TIP.2004.823819
% 
% Dependencies:
%   derivativesByFilter 
% %

% get some constants
dx = x(2)-x(1);
dy = y(2)-y(1);

% compute gradient
[dzdx, dzdy] = derivativesByFilters(zz, 'x', 'y', dx, dy, 'five');
