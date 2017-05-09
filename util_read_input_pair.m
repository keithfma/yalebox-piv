function [xx, yy, ini, ini_roi, fin, fin_roi] = util_read_input_pair(input_file, ini_step, gap)
% function [xw, yw, ini, ini_roi, fin, fin_roi] = util_read_input_pair(input_file, ini_step, gap)
%
% Reads pre-processed image pair from netCDF file produced by prep_series().
% Useful for parameter testing the piv() routine.
%
% Arguments:
%   input_file: String, file name of netCDF data file, as produced by
%       prep_series().
%   ini_step: Scalar, integer, step number of initial image frame in the image
%       pair to be returned.
%   gap: Scalar, steps between initial and final image in pair, defaults to 1 if
%       empty or not provided.
%
%   xx, yy = Vector, world coordinate vectors for image matrices
%   ini, fin = 2D matrix, intensity values for initial and final images 
%   ini_roi, fin_roi = 2D matrix, logical flags for inital and final images,
%       indicating if a pixel is sand (1) or not (0)
% %

%% check inputs

narginchk(2, 3);

if nargin < 3 || isempty(gap)
    gap = 1;
end

validateattributes(input_file, {'char'}, {'vector'});
validateattributes(ini_step, {'numeric'}, {'scalar', 'integer', '>=', 0});
validateattributes(gap, {'numeric'}, {'scalar', 'integer', 'positive'});

%% read data

% get image coordinate vectors
xx = double(ncread(input_file, 'x'));
yy = double(ncread(input_file, 'y'));

% get constant mask
mask_manual = logical(ncread(input_file, 'mask_manual'));

% find index of ini and fin
step = ncread(input_file, 'step');
ini_index = find(step == ini_step);
fin_index = ini_index + gap;

ini = double(ncread(input_file, 'img', [1, 1, ini_index], [inf, inf, 1]));
fin = double(ncread(input_file, 'img', [1, 1, fin_index], [inf, inf, 1]));

ini_mask = ncread(input_file, 'mask_auto', [1, 1, ini_index], [inf, inf, 1]); 
ini_roi = logical(ini_mask) & mask_manual;

fin_mask = ncread(input_file, 'mask_auto', [1, 1, fin_index], [inf, inf, 1]); 
fin_roi = logical(fin_mask) & mask_manual;