function [dzdx, dzdy] = post_gradient(x, y, zz, pad_method)
%
% Return numerical gradient using optimal 7-tap method from reference [1].
%
% Arguments:
%   x, y: Coordinate vectors for x- and y-directions
%   zz: Data matrix for which gradient is to be computed
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

% pad coordinate vectors
pad_coord = @(z, dz) [z(1)+dz*(-pad_width:-1)'; z(:); z(end)+dz*(1:pad_width)'];
x_p = pad_coord(x, x(2)-x(1));
y_p = pad_coord(y, y(2)-y(1));

% pad data matrix
if strcmp(pad_method, 'nearest')
    
    keyboard
    
else
    error('invalid padding method selected');
end
