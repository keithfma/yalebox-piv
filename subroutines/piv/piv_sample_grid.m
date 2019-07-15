function [r_grid, c_grid, pad_nr, pad_nc] = piv_sample_grid(sample_len, sample_spc, img_nr, img_nc)
%function [r_grid, c_grid, pad_nr, pad_nc] = piv_sample_grid(sample_len, sample_spc, img_nr, img_nc)
%
% Create sample grid and define image pad widths. The goal is a centered
% grid with a 2 extra observations at the edges, and all sample windows
% completely filled.
%
% Arguments:
%   sample_len = Vector, sample window size for each piv pass
%   sample_spc = Scalar, integer, grid spacing in pixels
%   img_nr, img_nc = Scalar, dimensions (number of rows and columns) of the
%       input images.
%
% Returns:
%   r_grid, c_grid = 2D matrix, integer, coordinate matrices, as constructed
%       by meshgrid, for the sample grid in the y (a.k.a. row) and x (a.k.a.
%       column) directions, in pixels
%   pad_nr, pad_nc = scalars, number of pixels to pad images (and thier
%       coordinates) in the row and column directions
% %

% validate inputs
validateattributes(sample_len, {'numeric'}, {'vector', 'integer'});
validateattributes(sample_spc, {'numeric'}, {'scalar', 'integer'});
validateattributes(img_nr, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(img_nc, {'numeric'}, {'scalar', 'integer', 'positive'});

% note: intentionally verbose to make this easier to understand

% number of "extra" samples to pad out in each direction
% note: need at least 1 for code to work, 3 is sufficient for 7-tap derivative 
num_extra = 3; 
assert(num_extra >= 1, 'constant num_extra breaks assumptions used below'); 

remainder = mod((img_nr - 1), sample_spc);
r_vec = (1 + remainder/2):sample_spc:img_nr;
r_vec = [...  % add extra samples
    (r_vec(1) - num_extra*sample_spc):sample_spc:r_vec(1)-sample_spc, ...
    r_vec, ...
    (r_vec(end) + sample_spc):sample_spc:(r_vec(end) + num_extra*sample_spc)] + num_extra*sample_spc;
r_vec = r_vec + sample_len(1)/2; % shift so first window is covered
pad_nr = num_extra*sample_spc + sample_len(1)/2;

remainder = mod((img_nc - 1), sample_spc);
c_vec = (1 + remainder/2):sample_spc:img_nc;
c_vec = [...  % add extra samples
    (c_vec(1) - num_extra*sample_spc):sample_spc:c_vec(1)-sample_spc, ...
    c_vec, ...
    (c_vec(end) + sample_spc):sample_spc:(c_vec(end) + num_extra*sample_spc)] + num_extra*sample_spc;
c_vec = c_vec + sample_len(1)/2; % shift so first window is covered
pad_nc = num_extra*sample_spc + sample_len(1)/2;

% shift grid to ensure sample windows span integer-pixel range
if is_even(sample_len)  % validated such that all even or all odd
    if is_whole(r_vec(1))
        r_vec = r_vec + 0.5;
    end
    if is_whole(c_vec(1))
        c_vec = c_vec + 0.5;
    end
elseif is_odd(sample_len)
    if ~is_whole(r_vec(1))
        r_vec = r_vec + 0.5;
    end
    if ~is_whole(c_vec(1))
        c_vec = c_vec + 0.5;
    end
else
    error(['sample_len must be either all-even or all-odd so that '
        'sample window is centered on the specified point']);
end

% bugfix: confirm that grids are spaced as expected, this may seem
%   obviously true given the code above, but it is critical and has been
%   broken in the past
assert(all(diff(r_vec) == sample_spc), 'row grid has incorrect spacing');
assert(all(diff(c_vec) == sample_spc), 'column grid has incorrect spacing');

% generate grid from vectors
[c_grid, r_grid] = meshgrid(c_vec, r_vec);

% % <DEBUG>: Check grid edges to confirm it is centered in the domain
% fprintf('Left edge [pixel]: %f, Right edge [pixel]: %f\n', ...
%     c_vec(1) - 1, img_nc - c_vec(end));
% fprintf('Left edge [world]: %f, Right edge [world]: %f\n', ...
%     x_sample_vector(1) - x_world(1), x_world(end) - x_sample_vector(end));
% fprintf('Bottom edge [pixel]: %f, Top edge [pixel]: %f\n', ...
%     r_vec(1) - 1, img_nr - r_vec(end));
% fprintf('Bottom edge [world]: %f, Top edge [world]: %f\n', ...
%     y_sample_vector(1) - y_world(1), y_world(end) - y_sample_vector(end));
% % </DEBUG>


function result = is_even(x)
result = all(mod(x, 2) == 0);


function result = is_odd(x)
result = all(mod(x, 2) == 1);


function result = is_whole(x)
result = mod(x, 1) == 0;
