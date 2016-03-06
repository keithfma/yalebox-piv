function [param] = default_movie_displ()
% function [param] = default_movie_displ()
%
% Generate default input struct for movie_displ() function. Typical usage
% would be to run this program to geneate a default parameter struct, then
% manually edit the members to customize.
% 
% See movie_displ() for information about each parameter.
% %

param.xunits = 'm';
param.yunits = 'm';
param.xlim = [-inf, inf];
param.ylim = [-inf, inf];
param.clim = [0.05, 0.95];
param.qsize = [20, 10]; 
param.qbnd = 0.05;
param.qscale = 01.;
param.tmp_dir = './tmp_movie_displacement';
param.tmp_file = fullfile(tmp_dir, 'tmp_%04i.png');
param.bbox = [];