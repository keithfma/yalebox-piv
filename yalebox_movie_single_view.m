function [] = yalebox_movie_single_view(input_file, output_file)
%
% input_file = String, path to input netCDF file as produced by
%   yalebox_prep_create_input().
%
% output_file = String path to output video file, without file extension
%
% Keith Ma, August 2015

% define parameters

% get netcdf ids
ncid = netcdf.open(input_file, 'NOWRITE');
x_id = netcdf.inqVarID(ncid, 'x');
y_id = netcdf.inqVarID(ncid, 'y');
step_id = netcdf.inqVarID(ncid, 'step');
intens_id = netcdf.inqVarID(ncid, 'intensity');

% get coordinate vectors
x = netcdf.getVar(ncid, x_id);
y = netcdf.getVar(ncid, y_id);
step = netcdf.getVar(ncid, step_id);
nx = numel(x);
ny = numel(y);
ns = numel(step);

% % prepare movie object: 16-bit grayscale Motion JPEG 2000, lossless compression
% movie_writer = VideoWriter(output_file, 'Archival');
% movie_writer.MJ2BitDepth = 16;
% movie_writer.FrameRate = 10; % frames/second
% movie_writer.open();

% loop: read data, annotate images, create frames
%for i = 1:ns
for i = 1    
    % read intensity image
    intens = netcdf.getVar(ncid, intens_id, [0, 0, i-1], [nx, ny, 1])';
    
    % add annotations using insertShape and similar
    
    % create frame    
end

% finalize
netcdf.close(ncid);
% movie_writer.close();

keyboard