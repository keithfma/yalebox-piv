function param = movie_displacement_default()
%
% Creates a struct containing all the parameters needed by movie_single(),
% populated with default values. Typical usage would be to call this function to
% get a properly formated parameter struct, edit the parameters as needed, then
% pass the edited struct to movie_single() to make the movie.
% 
% Arguments:
%   
%   param = Struct, see movie_displacement help for more information 
% %

% plot_displacement parameters (axis limits, quiver grid)
param.bbox = [];
param.xlim = [-inf, inf];
param.ylim = [-inf, inf];
param.ulim = [-inf, inf];
param.vlim = [-inf, inf];
param.mlim = [-inf, inf];
param.qsize = [20, 10];  
param.qbnd = 0.05;
param.qscale = 0.1;

% video
prm.frame_rate = 7;
prm.max_dim = [1920, 1080];
