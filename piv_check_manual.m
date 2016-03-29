function [] = piv_check_manual(image_file, displ_file, step, result_file)
%
% Estimate displacements manually by selecting matching points in an image pair.
% Writes results to a csv file with columns: row_ti, column_ti, row_tf,
% column_tf, row_tm, column_tm, u_chk_tm, v_chk_tm, u_piv_tm, v_piv_tm.
%
% Arguments:
% 
% image_file = String, path to netCDF file containing PIV input images, as
%   produced by prep_series.m
%
% displ_file = String, path to netCDF file containing PIV output data, as
%   produced by piv_series.m
%
% step = Scalar, step to be checked, this should be a valid time step in
%   displ_file, and so will probably at the midpoint between two images (e.g.
%   11.5)
%
% result_file = String, path to the output CSV file, the program will prompt
%   before overwriting
%
% Plan:
% 
% GUI with velocity magnitude, ini and fin in subplots, linked axes
% Button: add point, use impoint() to select a point in both images
% Button: Compute results: gather data from all the points, compute displacements and write the output csv
% Button: Restart from results file: read a results file and add impoint objects to restart
% %

%% check for sane inputs

%% create the GUI

%% Interactive point selection

%% Compute and write results

%% Restart from results file