function [] = post_series(piv_netcdf, post_netcdf, pro_bbox, retro_bbox)
% 
% Run post-processing analyses on PIV data, and saves the results and
% metadata in a netCDF file.
%
% Arguments:
%   piv_netcdf: string, file name of netCDF containing PIV results as
%       produced by piv_series.m
%   post_netcdf: string, file name of netCDF to create for post-processing
%       results data and metadata. Overwriting is not allowed to prevent
%       unhappy accidents
%   pro_bbox, retro_bbox: 4-element position vectors [xmin, ymin, width,
%       height] for the bounding boxes used to estimate pro- and retro-
%       plate displacements from PIV data. 
%
% Output netCDF is self-documenting
%
% % Keith Ma

% parse input arguments
validateattributes(piv_netcdf, {'char'}, {'nonempty'}, ...
    mfilename, 'piv_netcdf');
assert(exist(piv_netcdf, 'file') == 2, ...
    sprintf('input file %s does not exist', piv_netcdf));
validateattributes(post_netcdf, {'char'}, {'nonempty'}, ...
    mfilename, 'post_netcdf');
assert(exist(post_netcdf, 'file') ~= 2, ...
    sprintf('output file %s already exists', piv_netcdf));    
validateattributes(pro_bbox, {'numeric'}, {'vector', 'numel', 4}, ...
    mfilename, 'pro_bbox');
validateattributes(retro_bbox, {'numeric'}, {'vector', 'numel', 4}, ...
    mfilename, 'retro_bbox');

% read in PIV data from netCDF
%...stored
xx   = double(ncread(piv_netcdf, 'x')); 
yy   = double(ncread(piv_netcdf, 'y')); 
step = double(ncread(piv_netcdf, 'step'));
uu   = double(ncread(piv_netcdf, 'u'));
vv   = double(ncread(piv_netcdf, 'v'));
roi  = double(ncread(piv_netcdf, 'roi'));
%...derived
mm = sqrt(uu.^2+vv.^2);
num_steps = length(step);

% create output file

% run analyses for each timestep
uv_pro   = nan(num_steps, 2);
uv_retro = nan(num_steps, 2);

for ii = 1:num_steps
    us = uu(:,:,ii);
    vs = vv(:,:,ii);
    uv_pro(ii,:)   = post_displ_rect(xx, yy, us, vs, pro_bbox);
    uv_retro(ii,:) = post_displ_rect(xx, yy, us, vs, retro_bbox);
    
    % <DEBUG>
    fprintf('step %d: uv_pro: [%.4f, %.4f] m/step\n', ...
        ii, uv_pro(ii,:)*1000);
    fprintf('step %d: uv_retro: [%.4f, %.4f] m/step\n', ...
        ii, uv_retro(ii,:)*1000);
    %</DEBUG>
    
    % area budget
    % infinitesimal strain analysis
    % finite strain analysis
end

% finalize

% debug
keyboard