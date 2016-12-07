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

% define constants
pad_width = 3;

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

% compute spatial gradients
[dudx, dudy] = post_gradient(x, y, uu, pad_method);
