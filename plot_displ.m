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

keyboard