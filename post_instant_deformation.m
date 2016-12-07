function [] = post_instant_deformation(...
    x, y, uu, vv, roi, pad_method, verbose)
%
% Compute instantaneous/infinitesimal deformation parameters for input
% velocity fields.
%
% NOTE: Assumes x and y are regularly spaced (i.e. dx and dy are constant)
%
% Arguments:
%   x, y: Coordinate vectors for x- and y-directions
%   uu, vv: Velocity field matrices for x- and y-direction, units [m/step]
%   roi: Region-of-interest mask matrix (false for no data)
%   pad_method: (temporary) Select padding method, valid options are 
%       {'nearest'}
%   verbose: (optional) Flag indicating whether to print verbose messages 
%       (true) or not (false), default = false
%
% Modified from original version (deformation.m) by Mark Brandon. Changes
% include: ...
%
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
if nargin == 5
    verbose = false;
end
validateattributes(verbose, {'numeric', 'logical'}, {'scalar'}, ...
    mfilename, 'verbose');

if verbose
   fprintf('%s: start\n', mfilename); 
end

% TODO: Review the following calculations in detail

% calculate deformation-gradient tensor 
% NOTE: F = E + 1, where E is the displacement-gradient tensor
[F11, F12] = post_spatial_gradient(x, y, uu, pad_method);
[F21, F22] = post_spatial_gradient(x, y, vv, pad_method);
F11 = F11 + 1;
F22 = F22 + 1;

% calculate spin tensor (W) and stretch tensor (D) in its principal form
% NOTE: D1 and D2 are the principal extension and shortening rates,
%   respectively. theta contains the principal shortening rate directions
%   in radians
D1 = nan(size(uu));
D2 = nan(size(uu));
D2x = nan(size(uu));
D2y = nan(size(uu));
spin = nan(size(uu));
for k = 1:numel(uu)
    F = [F11(k), F12(k); F21(k),  F22(k)];
    if any(isnan(F(:)))
        % skip if F contains one or more nans
        continue;
    end
    % polar decomposition, F = VR, B = V^2 = FF', and R = (V^-1)F, where B
    %   is the Left Cauchy-Green tensor.
    B = F*F';
    % calculate eigen solution and sort with eigenvalues in descending order
    [T, lambda] = eig(B);
    [lambda, order] = sort(diag(lambda), 'descend');
    S = sqrt(lambda);
    D1(k) = log(S(1));
    D2(k) = log(S(2));
    T = T(:,order);
    % calculate unit vector for maximum shortening rate direction
    theta = acos(T(1,2));
    D2x(k) = cos(theta);
    D2y(k) = sin(theta);
    % calculate spin
    invV = T*diag(lambda.^(-1/2))*T';
    R = invV*F;
    spin(k) = atan2(R(2,1),R(1,1));
end

% calculate invariants, assuming plane strain (D3=0)
% Dt = sqrt(D1.^2 + D2.^2);
Dd = sqrt(((D1-D2).^2 + D1.^2 + D2.^2)/3);
Dv = sqrt(1/3)*(D1 + D2);

% calculate kinematic numbers
% Wk = 2*spin./(sqrt(2)*Dt);
% Ak = Dv./(sqrt(2)*Dt);
WkStar = 2*spin./(sqrt(2)*Dd);
AkStar = Dv./(sqrt(2)*Dd);

keyboard