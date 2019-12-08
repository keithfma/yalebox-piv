function [r_grid, c_grid] = piv_sample_grid(sample_len, sample_spc, x_vec, y_vec)
%function [r_grid, c_grid] = piv_sample_grid(sample_len, sample_spc, x_vec, y_vec)
%
% Create sample grid registered to the world coordinate origin
%
% Arguments:
%   sample_len = Vector, sample window size for each piv pass
%   sample_spc = Scalar, integer, grid spacing in pixels
%   x_vec, y_vec = x (columns) and y (rows) coordinate vectors for input images.
%
% Returns:
%   r_grid, c_grid = 2D matrix, integer, coordinate matrices, as constructed
%       by meshgrid, for the sample grid in the y (a.k.a. row) and x (a.k.a.
%       column) directions, in pixels
% %

% validate inputs
validateattributes(sample_len, {'numeric'}, {'vector', 'integer'});
validateattributes(sample_spc, {'numeric'}, {'scalar', 'integer'});
validateattributes(x_vec, {'numeric'}, {'vector'});
validateattributes(y_vec, {'numeric'}, {'vector'});

% note: intentionally verbose to make this easier to understand

% find the origin pixel
r_origin = find(abs(y_vec) < eps);
c_origin = find(abs(x_vec) < eps);
assert(numel(r_origin) == 1, 'Failed to find y-component of origin pixel');
assert(numel(c_origin) == 1, 'Failed to find x-component of origin pixel');

r_vec = [...
    fliplr((r_origin - sample_spc):-sample_spc:(sample_len(1)/2)), ...     % below origin
    r_origin, ...                                                          % origin
    (r_origin + sample_spc):sample_spc:(length(y_vec) - sample_len(1)/2)]; % above origin

assert(all(diff(r_vec) == sample_spc), 'sample grid spacing in y-direction is messed up');
assert(sum(r_vec == r_origin) == 1, 'sample grid does not contain origin in y-direction');
assert(all(r_vec >= 1) && all(r_vec <= (numel(y_vec) + 1)), 'sample grid out-of-bounds in y-direction');

c_vec = [...
    fliplr((c_origin - sample_spc):-sample_spc:(sample_len(1)/2)), ...     % below origin
    c_origin, ...                                                          % origin
    (c_origin + sample_spc):sample_spc:(length(x_vec) - sample_len(1)/2)]; % above origin

assert(all(diff(c_vec) == sample_spc), 'sample grid spacing in x-direction is messed up');
assert(sum(c_vec == c_origin) == 1, 'sample grid does not contain origin in x-direction');
assert(all(c_vec >= 1) && all(c_vec <= (numel(x_vec) + 1)), 'sample grid out-of-bounds in x-direction');

% generate grid from vectors
[c_grid, r_grid] = meshgrid(c_vec, r_vec);
