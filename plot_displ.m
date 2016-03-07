function plot_displ(piv_file, index, varargin)
% function plot_displ(piv_file, index, varargin)
%
% Plot displacement magnitude and direction. Generates 3 subplots:
%
%   - x-direction displacement 
%   - y-direction displacement 
%   - total displacement magnitude and direction
% 
% Arguments:
% Parameters (Name-Value pairs):
% %

%% parse input arguments

% positional arguments
validateattributes(piv_file, {'char'}, {'vector'});
assert(exist(piv_file, 'file')==2);
validateattributes(index, {'numeric'}, {'scalar', 'integer'});

% parameter name-value pairs
ip = inputParser();

ip.addParameter('coord_units', 'm', ...
    @(x) ismember(x, {'m', 'cm'})); 
ip.addParameter('displ_units', 'm/step', ...
    @(x) ismember(x, {'m/step', 'mm/step', '1'}));
ip.addParameter('norm_bbox', [], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 4}));
ip.addParameter('xlim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('ylim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('ulim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('vlim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('mlim', [-inf, inf], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
ip.addParameter('qsize', [20, 10], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2, 'integer'}));
ip.addParameter('qbnd', 0.05, ...
    @(x) validateattribute(x, {'numeric'}, {'scalar', '>=', 0, '<=', 0.5}));
ip.addParameter('qscale', 0.2, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));

ip.parse(varargin{:});
options = ip.Results;

%% prepare data

[step, xx, yy, uu, vv, mm] = util_read_piv_step(piv_file, index);

% convert coordinate units
switch options.coord_units
    case 'm'
        % default, do nothing
    case 'cm'
        xx = xx*100;
        yy = yy*100;
end

% convert displacement units
switch options.displ_units
    case 'm/step'
        % default, do nothing
    case 'mm/step'
        uu = uu*1000;
        vv = vv*1000;
        mm = mm*1000;
    case '1'
        [uu, vv, mm, options.norm_bbox] = ...
            util_normalize_displ(xx, yy, uu, vv, mm, options.norm_bbox);
end

%% plot



%% subroutines

keyboard