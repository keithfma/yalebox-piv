function [r_sample, c_sample] = piv_sample_grid(sample_len, sample_spc, img_nr, img_nc)
%function [r_sample, c_sample] = piv_sample_grid(sample_len, sample_spc, img_nr, img_nc)
%
% Create sample grid which is the largest centered grid that fits in the domain
%
% Arguments:
%   sample_len = Vector, sample window size for each piv pass
%   sample_spc = Scalar, integer, grid spacing in pixels
%   img_nr, img_nc = Scalar, dimensions (number of rows and columns) of the
%       input images.
%   r_sample, c_sample = 2D matrix, integer, coordinate matrices, as constructed
%       by meshgrid, for the sample grid in the y (a.k.a. row) and x (a.k.a.
%       column) directions, in pixels
% %

% TODO: validate inputs
% TODO: validate samp_len is all even or all odd, required so that sample
%   window is centered on the specified point

remainder = mod((img_nr - 1), sample_spc);
r_sample_vector = (1 + remainder/2):sample_spc:img_nr;

remainder = mod((img_nc - 1), sample_spc);
c_sample_vector = (1 + remainder/2):sample_spc:img_nc;

% shift grid to ensure sample windows span integer-pixel range
if is_even(sample_len(1))  % validated such that all even or all odd
    if is_whole(r_sample_vector(1))
        r_sample_vector = r_sample_vector + 0.5;
    end
    if is_whole(c_sample_vector(1))
        c_sample_vector = c_sample_vector + 0.5;
    end
else
    if ~is_whole(r_sample_vector(1))
        r_sample_vector = r_sample_vector + 0.5;
    end
    if ~is_whole(c_sample_vector(1))
        c_sample_vector = c_sample_vector + 0.5;
    end
end

% generate grid from vectors
[c_sample, r_sample] = meshgrid(c_sample_vector, r_sample_vector);

% % <DEBUG>: Check grid edges to confirm it is centered in the domain
% fprintf('Left edge [pixel]: %f, Right edge [pixel]: %f\n', ...
%     c_sample_vector(1) - 1, img_nc - c_sample_vector(end));
% fprintf('Left edge [world]: %f, Right edge [world]: %f\n', ...
%     x_sample_vector(1) - x_world(1), x_world(end) - x_sample_vector(end));
% fprintf('Bottom edge [pixel]: %f, Top edge [pixel]: %f\n', ...
%     r_sample_vector(1) - 1, img_nr - r_sample_vector(end));
% fprintf('Bottom edge [world]: %f, Top edge [world]: %f\n', ...
%     y_sample_vector(1) - y_world(1), y_world(end) - y_sample_vector(end));
% % </DEBUG>


function result = is_even(x)
result = mod(x, 2) == 0;


function result = is_whole(x)
result = mod(x, 1) == 0;
