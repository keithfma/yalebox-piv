% Skeleton for re-implementation of yalebox PIV analysis routine

% Arguments:
%
%   ini = 2D matrix, double, range 0 to 1, normalize grayscale image from
%       the start of the step to be analyzed.
%
%   fin = 2D matrix, double, range 0 to 1, normalize grayscale image from
%       the end of the step to be analyzed.
%
%   x = Vector, double, increasing, x-direction coordinate vector, length
%       must match the columns in ini and fin.
%
%   y = Vector, double, increasing, y-direction coordinate vector, length
%       must match the rows in ini and fin.
%
%   npass = Scalar, integer, number of PIV grid refinement passes
%
%   samplen = Vector, length === npass, double, approximate side length of
%       the sample window in world coordinates (same as x and y),
%       approximate-only because it must be rounded to the nearest odd
%       integer in pixels.
%
% References:


% debug parameters
data_dir = '/home/kfm/Documents/dissertation/yalebox-exp-fault/data/fault_ss_01/piv/';

% input arguments, hard-coded for debug
ini = flipud(double(imread([data_dir, 'img1.png']))/double(uint16(inf)));
fin = flipud(double(imread([data_dir, 'img2.png']))/double(uint16(inf)));
x = (1:size(ini,2))/1e3;
y = (1:size(ini,1))/1e3;
npass = 1;
sampwid = 0.1;

% validate input arguments
validateattributes(ini,...
    {'double'}, {'2d', 'real', 'nonnan', '>=', 0, '<=' 1}, ...
    mfilename, 'ini');
validateattributes(fin,...
    {'double'}, {'2d', 'real', 'nonnan', '>=', 0, '<=' 1}, ...
    mfilename, 'fin');
validateattributes(x, ...
    {'double'}, {'vector', 'real', 'nonnan', 'increasing'}, ...
    mfilename, 'x');
validateattributes(y, ...
    {'double'}, {'vector', 'real', 'nonnan', 'increasing'}, ...
    mfilename, 'y');
validateattributes(npass, ...
    {'numeric'}, {'scalar', 'integer', 'nonnegative'}, ...
    mfilename, 'npass');
validateattributes(sampwid, ...
    {'double'}, {'numel', npass, 'real', 'nonnan', 'positive'}, ...
    mfilename, 'sampwid');

% grid refinement loop (start)

% convert lengths from world to pixels

% compute sample window center points in pixels coords

% pad data and sample window center points (window+search_range+offsets, use padarray)

% (can I normalize away the penalty for flow out of the image?)


