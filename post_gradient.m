function [dzdx, dzdy] = post_gradient(x, y, zz, pad_method)
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
%
% Dependencies:
%
% %

% define parameters
pad_width = 3;

% check inputs
narginchk(4,4);
validateattributes(x, {'numeric'}, {'vector','real'}, mfilename, 'x');
validateattributes(y, {'numeric'}, {'vector','real'}, mfilename, 'y');
nx = length(x);
ny = length(y);
validateattributes(zz, {'numeric'}, {'size', [ny, nx]}, mfilename, 'zz');
validateattributes(pad_method, {'char'}, {'vector'}, mfilename, 'pad_method');

% pad coordinate vectors and data matrix
pad_coord = @(z, dz) [z(1)+dz*(-pad_width:-1)'; z(:); z(end)+dz*(1:pad_width)'];
x_p = pad_coord(x, x(2)-x(1));
y_p = pad_coord(y, y(2)-y(1));

[xx_p, yy_p] = meshgrid(x_p, y_p);
zz_p = padarray(zz, [pad_width, pad_width], NaN, 'both');

if strcmp(pad_method, 'nearest')
    roi = ~isnan(zz_p);
    si = scatteredInterpolant(xx_p(roi), yy_p(roi), zz_p(roi), ...
        'nearest', 'nearest');
    zz_p(~roi) = si(xx_p(~roi), yy_p(~roi));
else
    error('invalid padding method selected');
end

keyboard
