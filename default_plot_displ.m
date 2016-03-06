function [param] = default_plot_displ(is_norm)
% function [param] = default_plot_displ(is_norm)
%
% Generate default input struct for plot_displ() function. Typical usage would
% be to run this program to geneate a default parameter struct, then manually
% edit the members to customize the plot.
%
% Arguments:
%
%   norm= Scalar, logical, use value for 'normalized' plot (1), or don't (0)
%
% See plot_displ() for information about each parameter.
% %

param.xunits = 'm';
param.yunits = 'm';
if is_norm    
    param.uunits = '1';
    param.vunits = '1';
    param.munits = '1';
else
    param.uunits = 'm/step';
    param.vunits = 'm/step';
    param.munits = 'm/step';
end
param.xlim = [-inf, inf];
param.ylim = [-inf, inf];
param.ulim = [-inf, inf];
param.vlim = [-inf, inf];
param.mlim = [-inf, inf];
param.qsize = [20, 10]; 
param.qbnd = 0.05;
param.qscale = 01.;