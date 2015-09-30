function [xx, yy, uu, vv] = ...
    yalebox_piv_step(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0, verbose)
%
% function [xx, yy, uu, vv] = ...
%     yalebox_piv_step(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0. verbose)
%
% New implementation PIV analysis for Yalebox image data
%
% Arguments, input:
%
%   ini = 2D matrix, double, range 0 to 1, normalized grayscale image from
%       the start of the step to be analyzed.
%
%   fin = 2D matrix, double, range 0 to 1, normalized grayscale image from
%       the end of the step to be analyzed.
%
%   xx = Vector, double, increasing, x-direction coordinate vector, length
%       must match the columns in ini and fin.
%
%   yy = Vector, double, increasing, y-direction coordinate vector, length
%       must match the rows in ini and fin.
%
%   samplen = Scalar, integer, side length of the square sample window
%
%   sampspc = Scalar, integer, spacing between adjacent sample points in the
%       (square) sample grid
%
%   intrlen = Scalar, integer, side length of the square interrogation window
%
%   u0 = Scalar (PLANNED: or 2D matrix), integer, initial guess for x-direction
%       displacement, units are [pixels]
%
%   v0 = Scalar (PLANNED: or 2D matrix), integer, initial guess for y-direction
%       displacement, units are [pixels]
%
%   verbose = Scalar, integer, flag to enable (1) or diasable (0) verbose text
%       output messages
%
% Arguments, output:
%
%   xx, yy = Vector, double, coordinate vectors for the final output sample
%       grid, in world coordinate units
%
%   uu, vv = 2D matrix, double, computed displacement in the x- and y-directions
%       in world coordinate units
%
% References:
%
% [1] Raffel, M., Willert, C. E., Wereley, S. T., & Kompenhans, J. (2007).
%   Particle Image Velocimetry: A Practical Guide. BOOK, 1–448. 
%
% [2] Wereley, S. T. (2001). Adaptive Second-Order Accurate Particle Image
%   Velocimetry, 31
%

% parse inputs
check_input(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0, verbose);
if verbose
    print_sep('Input Arguments');
    print_input(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0);
end

% init sample grid 
[rr, cc] = sample_grid(sampspc, size(ini));
nr = length(rr);
nc = length(cc);

% init displacements
uu = u0*ones(nr, nc); % set to  initial guess
vv = v0*ones(nr, nc); 

% debug {
xx = 0;
yy = 0;
% } debug 

keyboard

end

%% computational subroutines

function [] = check_input(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0, verbose)
% Check for sane input argument properties, exit with error if they do not
% match expectations.
              
validateattributes(ini,...
    {'double'}, {'2d', 'real', 'nonnan', '>=', 0, '<=' 1}, ...
    mfilename, 'ini');
[nr, nc] = size(ini);

validateattributes(fin,...
    {'double'}, {'2d', 'real', 'nonnan', '>=', 0, '<=' 1, 'size', [nr, nc]}, ...
    mfilename, 'fin');

validateattributes(xx, ...
    {'double'}, {'vector', 'real', 'nonnan', 'numel', nc}, ...
    mfilename, 'xx');

validateattributes(yy, ...
    {'double'}, {'vector', 'real', 'nonnan', 'numel', nr}, ...
    mfilename, 'yy');

validateattributes(samplen, ...
    {'numeric'}, {'scalar', 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'samplen');

validateattributes(sampspc, ...
    {'numeric'}, {'scalar', 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'sampspc');

validateattributes(intrlen, ...
    {'numeric'}, {'scalar', 'integer', 'positive', 'nonnan', }, ...
    mfilename, 'intrlen');

validateattributes(u0, ...
    {'double'}, {'scalar'}, ...
    mfilename, 'u0');

validateattributes(v0, ...
    {'double'}, {'scalar'}, ...
    mfilename, 'u0');

validateattributes(verbose, {'numeric', 'logical'}, {'scalar', 'binary'}, ...
    mfilename, 'verbose');

end

function [rr, cc] = sample_grid(spc, sz)
% Create sample grid using equally spaced points at integer pixel
% coordinates, with the remainder split evenly between edges.
%
% Arguments:
%
%   spc = Scalar, integer, grid spacing in pixels
%
%   sz = Vector, length == 2, grid dimensions [rows, columns] for the
%       original input data matrices (e.g. ini and fin).
%
%   rr, cc = Vector, integer, coordinate vectors for the sample grid in the
%       y (a.k.a. row) and x (a.k.a. column) directions, in pixels

rrem = mod(sz(1)-1, spc);
rr0 = 1+floor(rrem/2);
rr = rr0:spc:sz(1);

crem = mod(sz(2)-1, spc);
cc0 = 1+floor(crem/2);
cc = cc0:spc:sz(2);

end

%% verbose subroutines --------------------------------------------------

function [] = print_sep(msg)
% Print a user-specified message and a separator line for verbose output
% messages

fprintf('----------\n%s\n', msg);

end

function print_input(ini, fin, xx, yy, samplen, sampspc, intrlen, u0, v0)
% Display values (or a summary of them) for the input arguments

fprintf('ini: size = [%i, %i], fraction data = %.2f%%\n',...
    size(ini, 1), size(ini, 2), sum(ini(:) ~= 0)/numel(ini)*100);

fprintf('fin: size = [%i, %i], fraction data = %.2f%%\n',...
    size(fin, 1), size(fin, 2), sum(fin(:) ~= 0)/numel(fin)*100);

fprintf('xx: length = %i, min = %.3f, max = %.3f, delta = %.3f\n', ...
    length(xx), min(xx), max(xx), xx(2)-xx(1));

fprintf('yy: length = %i, min = %.3f, max = %.3f, delta = %.3f\n', ...
    length(yy), min(yy), max(yy), yy(2)-yy(1));

fprintf('samplen: %s\n', sprintf('%i  ', samplen));

fprintf('sampspc: %s\n', sprintf('%i  ', sampspc));

fprintf('intrlen: %s\n', sprintf('%i  ', intrlen));

fprintf('u0: %i\n', u0);

fprintf('u0: %i\n', v0);

end