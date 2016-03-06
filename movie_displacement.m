function [] = movie_displacement(piv_file, movie_file, clim, xlim, ylim, qsize, qbnd, qscale, cleanup, show_frame)
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
%   cleanup = OPTIONAL, delete temporary image files (1) or don't (0)
%
%   show_frame = OPTIONAL, show a single frame at step == show_frame, do not
%       process other frames or make a movie, used for testing parameter values.
% %

% local constants
tmp_dir = './tmp_movie_displacement';

% parse inputs
% make_movie = nargin < 9 || isempty(show_frame);

% get color axis limits from vector quantiles
uu = ncread(piv_file, 'u');
vv = ncread(piv_file, 'v');
mm = sqrt(uu.^2 + vv.^2);
ulim = quantile(uu(:), clim);
vlim = quantile(vv(:), clim);
mlim = quantile(mm(:), clim);
num_steps = size(uu, 3);

% init loop
clear uu vv mm
bbox = [];
mkdir(tmp_dir);


% loop over all timesteps
for ii = 1:100; %num_steps
    
    % plot frame
    bbox = plot_displacement_norm(piv_file, ii, bbox, xlim, ylim, ulim, vlim, mlim, ...
        qsize, qbnd, qscale);
    
    % first time: get parameters needed to convert figure to FHD image (1920x1080)
    if ii == 1
        % get magnification factor from simple test image
        img_size = size(export_fig('-dpng'));
        img_size = img_size(1:2); % rows, cols
        magnify = sprintf('-m%.10f', min([1080, 1920]./img_size));       
        
        % get padding from magnified test image
        img_size = size(export_fig('-dpng', magnify));
        img_size = img_size(1:2);
        vpad = [floor((1080-img_size(1))/2), ceil((1080-img_size(1))/2)]; 
        hpad = [floor((1920-img_size(2))/2), ceil((1920-img_size(2))/2)];         
    end
    
    % convert figure to FHD image (1920x1080)
    img = export_fig('-dpng', magnify);
    img = padarray(img, [vpad(1), hpad(1), 0], 1, 'pre');
    img = padarray(img, [vpad(2), hpad(2), 0], 1, 'post');
    imwrite(img, sprintf('tmp_%04i.png', ii));

    close(gcf);
end
