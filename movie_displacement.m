function [] = movie_displacement(piv_file, movie_file, clim, xlim, ylim, qsize, qbnd, qscale, show_frame)
% Generate movie from PIV displacement results. Calls plot_displacement_norm()
% to create frames for each time step. 
%
% Arguments:
%
%   piv_file = String, path to netCDF containing PIV results produced by piv.m
%
%   movie_file = String path to output video file, without file extension
%
%   clim =  Vector, length ==2, color axis limits as quantiles of the whole data
%       set. Sets ulim, vlim, mlim parameters to plot_displacement
%   
%   xlim, ylim = 
%
%   qsize, qbnd, qscale = quiver plot parameters, see plot_displacement help for
%       details
%
%   show_frame = OPTIONAL, show a single frame at step == show_frame, do not
%       process other frames or make a movie, used for testing parameter values.
% %

% parse inputs
make_movie = nargin < 7 || isempty(show_frame);

% get color axis limits from vector quantiles
uu = ncread(piv_file, 'u');
vv = ncread(piv_file, 'v');
mm = sqrt(uu.^2 + vv.^2);
ulim = quantile(uu(:), clim);
vlim = quantile(vv(:), clim);
mlim = quantile(mm(:), clim);
num_steps = size(uu, 3);
clear uu vv mm

% loop over all timesteps
bbox = [];
for ii = 1:num_steps
    
    % plot frame
    bbox = plot_displacement_norm(piv_file, ii, bbox, xlim, ylim, ulim, vlim, mlim, ...
        qsize, qbnd, qscale);
    
    pause(1)
    close(gcf);

end
