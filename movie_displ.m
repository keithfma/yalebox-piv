function opt = movie_displ(piv_file, movie_file, varargin)
% Generate movie from PIV displacement results. 
%
% Arguments:
%
%   piv_file = String, path to netCDF containing PIV results produced by piv.m
%
%   movie_file = String path to output video file, without file extension.
%       Movies are always created in the mp4 format with h264 encoding.
%
% Parameters:
%
%   'show_frame' = show a single frame at index == show_frame, do not
%       process other frames or make a movie. Used for testing parameter values.
%
%   'coord_units' = String, name of units for coordinate axes, coordinate values
%       will be rescaled accordingly, and labels will reflect these units. Valid
%       options are: 'm', 'cm'. Default = 'cm'
%
%   'norm_bbox' = Vector, length==4, bounding box for the data region used in
%       computing the displacement unit normalization factor. This is only used
%       if displ_units == '1'. If empty, the normalization function will prompt
%       the user to select a box interactively the first time. Default = [].
%
%   'xlim', 'ylim' = Vector, length==2, [minimum, maximum] values for the x-axis
%       and y-axis. Values beyond the range of the data will be truncated to the
%       data limits. Thus, to span the data, one could use [-inf, inf]. Default
%       = [-inf, inf].
%
%   'clim' =  Vector, length==2, color axis limits as quantiles of the whole
%       data set. Sets ulim, vlim, mlim parameters to plot_displ()
%
%   'qsize' = Vector, length==2, size of the grid for vector direction overlay
%       (quiver) in [rows, cols]
%
%   'qbnd' = Scalar, boundary margin for vector direction overlay, as fraction
%       of the axis ranges.
%
%   'qscale' = Scalar, range [0,1], length of vector lines in vector direction
%       overlay 
%
%   'font_size_title', 'font_size_tick', 'font_size_axis' = Scalar,
%       integer, font size for titles, tick labels, and axis labels,
%       repectively. Defaults are 14, 12, 12.
%
%   'tmp_dir' = String, directory where frame images should be saved. Default =
%       './tmp_movie_displacement'
%
%   'tmp_file' = String, fprintf-style formatting string defining format of
%       frame image files, must contain one and only one integer variable.
%       Also used by ffmpeg to identify input images. Default = 'tmp_%04d.png'
%
%   'fps' = 
% %

%% parse input arguments

% positional arguments
validateattributes(piv_file, {'char'}, {'vector'});
assert(exist(piv_file, 'file') == 2);
validateattributes(movie_file, {'char'}, {'vector'});
assert(exist(movie_file, 'file') ~= 2); % do not overwrite

% parameter name-value pairs
ip = inputParser();

ip.addParameter('show_frame', 0, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
ip.addParameter('coord_units', 'cm', ...
    @(x) ismember(x, {'m', 'cm'})); 
ip.addParameter('norm_bbox', [], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 4}));
ip.addParameter('xlim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('ylim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('clim', [0.05, 0.95], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2, '>=', 0, '<=' 1}));
ip.addParameter('qsize', [20, 10], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2, 'integer'}));
ip.addParameter('qbnd', 0.05, ...
    @(x) validateattribute(x, {'numeric'}, {'scalar', '>=', 0, '<=', 0.5}));
ip.addParameter('qscale', 0.2, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
ip.addParameter('font_size_title', 14, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
ip.addParameter('font_size_tick', 12, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
ip.addParameter('font_size_axis', 12, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
ip.addParameter('tmp_dir', './tmp_movie_displacement', ...
    @(x) validateattributes(x, {'char'}, {'vector'}));
ip.addParameter('tmp_file', 'tmp_%04d.png', ...
    @(x) validateattributes(x, {'char'}, {'vector'}));
ip.addParameter('fps', 10, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));

ip.parse(varargin{:});
opt = ip.Results;

%% get color axis limits from global quantiles of normalized data

% read all data to memory 
xx = ncread(piv_file, 'x'); 
yy = ncread(piv_file, 'y'); 
uu = permute( ncread(piv_file, 'u'), [2, 1, 3] );
vv = permute( ncread(piv_file, 'v'), [2, 1, 3] );
mm = sqrt(uu.^2 + vv.^2);
num_steps = size(uu, 3);

% normalize each step
for ii = 1:num_steps    
    [~, ~, uu(:,:,ii), vv(:,:,ii), mm(:,:,ii), opt.norm_bbox] = ...
        util_convert_displ_units(xx, yy, uu(:,:,ii), vv(:,:,ii), mm(:,:,ii), ...
            opt.coord_units, '1', opt.norm_bbox);
end

% compute limits from quantiles
ulim = quantile(uu(:), opt.clim);
vlim = quantile(vv(:), opt.clim);
mlim = quantile(mm(:), opt.clim);

% clean up section
clear xx yy uu vv mm

%% generate frame(s)

mkdir(opt.tmp_dir);
have_size = 0;

for ii = 333:num_steps
    
    % test case: skip all but specified frame
    if opt.show_frame ~= 0 && opt.show_frame ~= ii
        continue
    end
    
    % create plot
    plot_displ(piv_file, ii, ...
        'coord_units', opt.coord_units, ...
        'displ_units', '1', ...
        'norm_bbox', opt.norm_bbox, ...
        'xlim', opt.xlim, ...
        'ylim', opt.ylim, ...
        'ulim', ulim, ...
        'vlim', vlim, ...
        'mlim', mlim, ...
        'qsize', opt.qsize, ...
        'qbnd', opt.qbnd, ...
        'qscale', opt.qscale, ...
        'font_size_title', opt.font_size_title, ...
        'font_size_tick', opt.font_size_tick, ...
        'font_size_axis', opt.font_size_axis);
    
    % first time: get parameters needed to convert figure to FHD image (1920x1080)
    if have_size == 0
        
        % get magnification factor from simple test image
        img_size = size(export_fig('-dpng'));
        img_size = img_size(1:2); % rows, cols
        fact = min([1080, 1920]./img_size);
        fact = 0.99*fact; % avoid overshooting size by reducing then and pading
        magnify = sprintf('-m%.10f', fact);  
        
        % get padding from magnified test image
        img_size = size(export_fig('-dpng', magnify));
        img_size = img_size(1:2);
        vpad = [floor((1080-img_size(1))/2), ceil((1080-img_size(1))/2)]; 
        hpad = [floor((1920-img_size(2))/2), ceil((1920-img_size(2))/2)]; 
        
        have_size = 1;
    end
    
    % resize to FHD (1920x1080), and write to image file
    img = export_fig('-dpng', magnify);
    img = padarray(img, [vpad(1), hpad(1), 0], 1, 'pre');
    img = padarray(img, [vpad(2), hpad(2), 0], 1, 'post');
    imwrite(img, fullfile(opt.tmp_dir, sprintf(opt.tmp_file, ii)));

    % test case: leave figure open
    if opt.show_frame == 0
        close(gcf);
    end
    
end

%% create movie from frame images

% external call to ffmpeg, command follows the guide at:
%   https://trac.ffmpeg.org/wiki/Create%20a%20video%20slideshow%20from%20images

cmd = sprintf('ffmpeg -framerate %d -i %s -c:v libx264 -pix_fmt yuv420p %s.mp4', ...
    opt.fps, fullfile(opt.tmp_dir, opt.tmp_file), movie_file);
system(cmd);

    
