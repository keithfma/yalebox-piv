function [] = movie_displ(piv_file, movie_file, param, show_frame)
% Generate movie from PIV displacement results. 
%
% Arguments:
%
%   piv_file = String, path to netCDF containing PIV results produced by piv.m
%
%   movie_file = String path to output video file, without file extension
%
%   caxis_quantile =  Vector, length==2, color axis limits as quantiles of the
%       whole data set. Sets ulim, vlim, mlim parameters to plot_displacement
%
%   plot_param = Struct. Formatting parameters for the plot_displ() function.
%       Note that the .ulim, .vlim, .mlim members are ignored. These values are
%       set adaptively based on the caxis_quantile argument.
%
%   cleanup = OPTIONAL, delete temporary image files (1) or don't (0)
%
%   show_frame = OPTIONAL, show a single frame at step == show_frame, do not
%       process other frames or make a movie, used for testing parameter values.
% %

% 

% get color axis limits from global quantiles of normalized data
%...read all data to memory 
uu = ncread(piv_file, 'u');
vv = ncread(piv_file, 'v');
mm = sqrt(uu.^2 + vv.^2);
%...normalize each step
num_steps = size(uu, 3);
for ii = 1:num_steps
    [uu(:,:,ii), vv(:,:,ii), mm(:,:,ii), param.bbox] = ...
        util_normalize_displ(xx, yy, uu(:,:,ii), vv(:,:,ii), mm(:,:,ii), param.bbox);     
end
%...compute limits from quantiles
plot_param.ulim = quantile(uu(:), param.clim);
plot_param.vlim = quantile(vv(:), param.clim);
plot_param.mlim = quantile(mm(:), param.clim);
clear uu vv mm

% movie param -> plot param

% loop over all timesteps
mkdir(tmp_dir);
for ii = 1:num_steps
    
    % plot frame
    plot_displ(piv_file, ii, bbox, xlim, ylim, ulim, vlim, mlim, ...
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
    imwrite(img, sprintf(tmp_file, ii));

    close(gcf);
end

if cleanup
    rmdir(tmp_dir, 's');
end

if 