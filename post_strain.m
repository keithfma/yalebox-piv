function [L, F] = post_strain(x, y, uu, vv, roi, pad_method)
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
%   L: Displacement-gradient tensor field, 3D matrix, last dimension stores
%       tensor components column-wise (L11, L21, L12, L22)
%   F: Deformation-gradient tensor, 3D matrix, stored same as L
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

% calculate displacement-gradient tensor 
[L11, L12] = spatial_gradient(x, y, uu, pad_method);
[L21, L22] = spatial_gradient(x, y, vv, pad_method);

% calculate deformation-gradient tensor 
% NOTE: F = L + I, where L is the displacement-gradient tensor 
F11 = L11 + 1;
F12 = L12;
F21 = L21;
F22 = L22 + 1;

% TODO: include additional strain parameters

% prepare outputs
L = cat(3, L11, L21, L12, L22); % tensor store comlum-wise in 3rd dim
F = cat(3, F11, F21, F12, F22);


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