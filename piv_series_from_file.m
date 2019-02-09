function [] = prep_series_from_file(result_file, param_file)
% function [] = prep_series_from_file(result_file, param_file)
% 
% Run PIV analysis for a given image series. Simple wrapper around
% piv_series() that reads the (many) input parameters from the standard
% JSON file rather than as arguments.
%
% Arguments:
%   result_file = String, filename of the MAT input file to be created. 
%   param_file: JSON file containing all input parameters
% %

% load dependencies
update_path('jsonlab');

% load parameters struct
param = loadjson(param_file, 'SimplifyCell', 1);


