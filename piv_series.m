function [] = piv_series(output_file, input_file, samplen, sampspc, intrlen, ...
                  npass, valid_max, valid_eps)
% 
% Run PIV analysis for a given input series. Input is expected to be a netCDF
% file as created by prep_series(). Results are saved in a new netCDF file which
% includes all relevant metadata.
%
% output_file = String, name of the output netCDF file containing PIV results
%
% input_file = String, name of the input netCDF file containing pre-processed
%   image data
%
% samplen, sampspc, intrlen, npass, valid_max, valid_eps = Select input
%   variables for PIV routine piv(), other inputs are contained in the input
%   file.
%
% %