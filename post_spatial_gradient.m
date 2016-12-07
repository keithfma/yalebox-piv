function [dzdx, dzdy] = post_spatial_gradient(x, y, zz, pad_method)
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

% get some constants
dx = x(2)-x(1);
dy = y(2)-y(1);
roi = ~isnan(zz);

% pad coordinate vectors and data matrix
pad_coord = @(z, dz) [z(1)+dz*(-pad_width:-1)'; z(:); z(end)+dz*(1:pad_width)'];
x_p = pad_coord(x, dx);
y_p = pad_coord(y, dy);
[xx_p, yy_p] = meshgrid(x_p, y_p);
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