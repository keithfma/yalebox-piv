function [] = movie_displacement(input_file, output_file, param, show_frame)
%
% input_file = String, path to netCDF containing PIV results produced by piv.m
%
% output_file = String path to output video file, without file extension
%
% param = Struct, movie parameters, as produced by movie_displacement_default.m
%   .xlim   -|
%   .ylim    |
%   .ulim    |
%   .vlim    |
%   .mlim    |-> See plot_displacement help for details
%   .qsize   |
%   .qbnd    |
%   .qscale -|
%   .frame_rate = Video frame rate in FPS
%   .max_dim = Maximum video dimensionsions in pixels
%
% show_frame = OPTIONAL, show a single frame at step == show_frame, do not
%   process other frames or make a movie, used for testing parameter
%   values.
% %

% parse inputs
make_movie = nargin < 4 || isempty(show_frame);