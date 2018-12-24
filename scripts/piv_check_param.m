% Test PIV parameters
%
% Expected variables:
%   PREP_PARAM_FILE: path to parameters definition file. Use
%       templates/piv.json as a starting point, and populate the variables
%       therein to suit your experiment.
%   IMAGE_FILE: path to pre-processed image file netCDF, as produced by
%       prep_series.m
%
% This script is designed to be run cell-by-cell. Each cell runs one step
% of the image preparation process. Run a cell, inspect the results, and
% edit the parameter file until you are satisfied with the results. Note
% that some cells depend on the results of previous cells.
% %

addpath(fullfile(pwd, '..'));
update_path('piv');


 