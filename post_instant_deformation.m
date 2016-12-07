function [] = post_instant_deformation(x, y, uu, vv, roi, verbose)
%
% Compute instantaneous/infinitesimal deformation parameters for input
% velocity fields.
%
% Arguments:
%   x, y: Coordinate vectors for x- and y-directions
%   uu, vv: Velocity field matrices for x- and y-direction, units [m/step]
%   roi: Region-of-interest mask matrix (false for no data)
%   verbose: (optional) Flag indicating whether to print verbose messages 
%       (true) or not (false), default = false
%
% Modified from original version (deformation.m) by Mark Brandon. Changes
% include: ...
%
% %

% check inputs
narginchk(5,6);
validateattributes(x, {'numeric'}, {'vector','real'}, mfilename, 'x');
validateattributes(y, {'numeric'}, {'vector','real'}, mfilename, 'y');
nx = length(x);
ny = length(y);
validateattributes(uu, {'numeric'}, {'size', [ny, nx]}, mfilename, 'uu');
validateattributes(vv, {'numeric'}, {'size', [ny, nx]}, mfilename, 'vv');
validateattributes(roi,  {'logical'}, {'size', [ny, nx]}, mfilename, 'roi');
if nargin == 5
    verbose = false;
end
validateattributes(verbose, {'numeric', 'logical'}, {'scalar'}, ...
    mfilename, 'verbose');

if verbose
   fprintf('%s: start\n', mfilename); 
end