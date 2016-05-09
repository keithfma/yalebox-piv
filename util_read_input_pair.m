function [xx, yy, ini, ini_roi, fin, fin_roi] = util_read_input_pair(input_file, index)
%
% Reads pre-processed image pair from netCDF file produced by
% prep_series(). Useful for parameter testing the piv() routine.
%
% Arguments:
%
% input_file = String, file name of netCDF data file, as produced by
%   prep_series().
%
% index = Scalar, integer, index of initial image frame in the image
%   pair to be returned.
%
% xx, yy = Vector, world coordinate vectors for image matrices
%
% ini, fin = 2D matrix, intensity values for initial and final images 
%
% ini_roi, fin_roi = 2D matrix, logical flags for inital and final images,
%   indicating if a pixel is sand (1) or not (0)
% %

%% check inputs
validateattributes(input_file, {'char'}, {'vector'});

assert(exist(input_file, 'file') == 2);
info = ncinfo(input_file);
step_idx = find(strcmp({info.Variables(:).Name}, 'step'));
assert(~isempty(step_idx));
num_steps = info.Variables(step_idx).Size;

validateattributes(index, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', num_steps-1});

%% read data
xx = double(ncread(input_file, 'x'));
yy = double(ncread(input_file, 'y'));

roi_const = ncread(input_file, 'mask_manual', [1, 1], [inf, inf])';

ini = double( ncread(input_file, 'intensity', ...
    [1, 1, index], [inf, inf, 1])' );
fin = double( ncread(input_file, 'intensity', ...
    [1, 1, index+1], [inf, inf, 1])' );

ini_roi = ncread(input_file, 'mask_auto', ...
    [1, 1, index], [inf, inf, 1])' & roi_const;
fin_roi = ncread(input_file, 'mask_auto', ...
    [1, 1, index+1], [inf, inf, 1])' & roi_const;

%% done