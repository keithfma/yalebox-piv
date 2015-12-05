function [uu1, vv1] = yalebox_piv_interp2d(r0, c0, uu0, vv0, r1, c1, method)
% function [uu1, vv1] = yalebox_piv_interp2d(r0, c0, uu0, vv0, r1, c1, method)
%
% Interpolate a vector field from one gridded coordinate system to another. Uses
% the griddedInterpolate class, which supports extrapolation for all available
% interpolation methods and allows for efficient reuse of the intepolant. 
%
% Arguments:
%
%   r0, c0 = Vector, double, row and column coordinate vectors for the input
%     vector field in uu0, vv0
%
%   uu0, vv0 - 2d matrix, double, input vector field
%
%   r1, c1 = Vector, double, row and column coordinate vectors for the output
%     vector field returned in uu1, vv1
%
%   method = (Optional) string selecting intepolation method, must be supported
%     by the MATLAB griddedInterpolant class (R2015a includes: linear, nearest,
%     cubic, spline), default value is 'spline'.
%
%   uu1, vv1 = 2d matrix, double, output (interpolated) vector field
%
% %

% set defaults
if nargin < 7
  method = 'spline';
end

% get coordinate matrices
[cc0, rr0] = meshgrid(c0, r0);
[cc1, rr1] = meshgrid(c1, r1);

% create and query interpolant 
interpolant = griddedInterpolant(rr0, cc0, uu0, method, method);
uu1 = interpolant(rr1, cc1);

interpolant.Values = vv0;
vv1 = interpolant(rr1, cc1);
