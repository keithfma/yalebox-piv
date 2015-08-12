function [] = yalebox_movie_single(input_file, output_file)
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
max_mov_dim = [1920, 1080];

tri_color = 'red';
tri_opacity = 0.7;

title_font_size = 72;
title_font_color = 'red';
title_box_color = 'white';
title_box_opacity = 1;

scale_upleft = [0.1, 0.8]; % normalized
scale_width = 0.1; % m
scale_height = 0.05; % m


function aimg = yalebox_movie_aux_scalebar(img, x, y, ulc, wid, ht, bclr, fclr, fsize)

% get netcdf ids
ncid = netcdf.open(input_file, 'NOWRITE');
intens_id = netcdf.inqVarID(ncid, 'intensity');

% get coordinate vectors
x = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'x'));
y = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'y'));
step = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'step'));

% % prepare movie object: 16-bit grayscale Motion JPEG 2000, lossless compression
% movie_writer = VideoWriter(output_file, 'Archival');
% movie_writer.MJ2BitDepth = 16;
% movie_writer.FrameRate = 10; % frames/second
% movie_writer.open();

% loop: read data, annotate images, create frames
%for i = 1:numel(step)
for i = 1    

    % read intensity image
    intens = netcdf.getVar(ncid, intens_id, [0, 0, i-1], [numel(x), numel(y), 1])';
    
    % resize and reshape image as needed
    [frame, xf, yf] = yalebox_movie_aux_resize(intens, x, y, max_mov_dim);
    [frame, yf] = yalebox_movie_aux_flip(frame, yf);
    
    % add annotations (triangles, scale, title, counter)
    frame = yalebox_movie_aux_spoint(frame, xf, yf, tri_tip, tri_len, ...
        tri_color, tri_opacity);
    frame = yalebox_movie_aux_title(frame, title_str, title_font_size, ...
        title_font_color, title_box_color, title_box_opacity);

    % add frame to movie object   
    
end

% finalize
netcdf.close(ncid);
% movie_writer.close();

keyboard
