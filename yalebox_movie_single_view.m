function [] = yalebox_movie_single_view(input_file, output_file)
%
% input_file = String, path to input netCDF file as produced by
%   yalebox_prep_create_input().
%
% output_file = String path to output video file, without file extension
%
% Keith Ma, August 2015

% debug arguments
tri_tip = [0, 0; -0.15, 0]; % 1 triangle per row, position in m
tri_len = 0.01; % m
title_str = 'Party Time';

% internal parameters
tri_color = 'red';
tri_opacity = 0.7;
title_font_size = 72;
title_font_color = 'red';
title_box_color = 'white';
title_box_opacity = 0;

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
nx = numel(x); dx = abs(x(1)-x(2));
ny = numel(y); dy = abs(y(1)-y(2));
ns = numel(step);

% make sure image top is at row 1
if y(1) < y(end)
    flip = true;
    y = y(end:-1:1);
else 
    flip = false;
end

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
    
    % make sure image top is at row 1
    if flip; intens = flipud(intens); end
    
    % add annotations (triangles, scale, title, counter)
    frame = yalebox_movie_add_spt(intens, x, y, tri_tip, tri_len, ...
        tri_color, tri_opacity);
    frame = yalebox_movie_add_title(frame, title_str, title_font_size, ...
        title_font_color, title_box_color, title_box_opacity);

    % add frame to movie object   
    
end

% finalize
netcdf.close(ncid);
% movie_writer.close();

keyboard
