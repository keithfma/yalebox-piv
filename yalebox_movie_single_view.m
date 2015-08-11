function [] = yalebox_movie_single_view(input_file, output_file)
%
% input_file = String, path to input netCDF file as produced by
%   yalebox_prep_create_input().
%
% output_file = String path to output video file, without file extension
%
% Keith Ma, August 2015

% define parameters

% open netcdf, get coordinate vectors
ncid = netcdf.open(input_nc, 'NOWRITE');
x = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'x'));
y = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'y'));
step = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'step'));

% prepare movie object: 16-bit grayscale Motion JPEG 2000, lossless compression
movie_writer = VideoWriter(output_file, 'Archival');
movie_writer.MJ2BitDepth = 16;
movie_writer.FrameRate = 10; % frames/second
movie_writer.open();

% loop: read data, annotate images, create frames
nframes = numel(step);
nframes = 1;
for i = 1:nframes
    % read intensity image    
    % add annotations using insertShape and similar
    % create frame    
end

% finalize movie
movie_writer.close();

keyboard