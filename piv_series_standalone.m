function [] = piv_series_standalone(param_file)
% function [] = piv_series_standalone(param_file)
%
% Wrapper function to run piv_series() as a compiled standalone executable.
% Includes special configuration needed to run on the BU SCC cluster batch
% system.
%
% Arguments:
%   param_file: String, .MAT file containing all piv_series() input arguments,
%       typically produced using a customized version of
%       piv_get_param_template.m
% % 

% load parameters from file
args = load(param_file);
print(args);
