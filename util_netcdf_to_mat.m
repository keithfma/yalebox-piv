function [] = util_netcdf_to_mat(file_netcdf, file_mat, clobber)
% function [] = util_netcdf_to_mat(file_netcdf, file_mat, clobber)
% 
% Convert data and attributes from netCDF to MATLAB mat file format. Each netCDF
% variable is converted to a struct in the mat file, with members that hold the
% data and attributes for that variable. Global attributes are stored in a
% variable named global_attr.
%
% Arguments:
% 
%   file_netcdf = String, name of input netCDF file, required
%
%   file_mat = String, name of output mat file, required
% 
%   clobber = Scalar flag, overwrite output file (1) or don't (0), default = 0
%
% %

% set defaults 
if nargin < 3
    clobber = 0;
end

% sanity check
validateattributes(file_netcdf, {'char'}, {'vector'}, mfilename, 'file_netcdf');
validateattributes(file_mat, {'char'}, {'vector'}, mfilename, 'file_mat');
validateattributes(clobber, {'numeric', 'logical'}, {'scalar', 'binary'}, mfilename, 'clobber');

if exist(file_netcdf, 'file') ~= 2
    error('Specified input file does not exist');
end

if clobber == 0 && exist(file_mat, 'file') == 2
    error('Specified output file exists, and overwriting is disabled');
end

% init
info = ncinfo(file_netcdf);
out = struct(); % holds all data and attributes to be written to mat file 

% copy global attributes to struct
out.global_attr = struct();
for ii = 1:length(info.Attributes)
    name = strrep(info.Attributes(ii).Name, ' ', '_');
    val = info.Attributes(ii).Value;
    out.global_attr.(name) = val;
end

% copy netCDF variables and thier attributes to structs
for ii = 1:length(info.Variables)
    
    % create struct for variable
    var_name = strrep(info.Variables(ii).Name, ' ', '_');
    out.(var_name) = struct();
    
    % get dimension names as csv string
    dims = [];
    for jj = 1:length(info.Variables(ii).Dimensions)
        name = strrep(info.Variables(ii).Dimensions(jj).Name, ' ', '_');
        dims = [dims ', ' name]; %#ok!
    end
    out.(var_name).dimensions = dims(3:end);
    
    % get attributes
    for jj = 1:length(info.Variables(ii).Attributes)
        name = strrep(info.Variables(ii).Attributes(jj).Name, ' ', '_');
        val = info.Variables(ii).Attributes(jj).Value;
        out.(var_name).(name) = val;
    end
    
    % get data 
    out.(var_name).value = ncread(file_netcdf, info.Variables(ii).Name);
    
end

% save results as a mat file
save(file_mat, '-struct', 'out');