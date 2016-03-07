function opt = movie_displ(piv_file, movie_file, show_frame, varargin)
% Generate movie from PIV displacement results. 
%
% Arguments:
%
%   piv_file = String, path to netCDF containing PIV results produced by piv.m
%
%   movie_file = String path to output video file, without file extension
%
%   show_frame = OPTIONAL, show a single frame at index == show_frame, do not
%       process other frames or make a movie. Used for testing parameter values.
%
% Parameters:
%   'coord_units' = String, name of units for coordinate axes, coordinate values
%       will be rescaled accordingly, and labels will reflect these units. Valid
%       options are: 'm', 'cm'. Default = 'cm'
%
%   'displ_units' = String, name of units for displacement vector components and
%       magnitude, values will be rescaled accordingly, and labels will reflect
%       these units. Valid options are: 'm', 'mm', '1'. Default is 'mm'.
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
% %

%% parse input arguments

% positional arguments
validateattributes(piv_file, {'char'}, {'vector'});
assert(exist(piv_file, 'file') == 2);
validateattributes(movie_file, {'char'}, {'vector'});
assert(exist(movie_file, 'file') ~= 2); % do not overwrite

% parameter name-value pairs
ip = inputParser();

ip.addParameter('coord_units', 'cm', ...
    @(x) ismember(x, {'m', 'cm'})); 
ip.addParameter('displ_units', 'mm/step', ...
    @(x) ismember(x, {'m/step', 'mm/step', '1'}));
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

ip.parse(varargin{:});
opt = ip.Results;

keyboard
