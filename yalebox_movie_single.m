function [] = yalebox_movie_single(input_file, output_file, show_frame)
%
% input_file = String, path to input netCDF file as produced by
%   yalebox_prep_create_input().
%
% output_file = String path to output video file, without file extension
%
% show_frame = OPTIONAL, show a single frame at step == show_frame, do not
%   process other frames or make a movie, used for testing parameter
%   values.
%
% Keith Ma, August 2015

% debug arguments
tri_tip = [0, 0; -0.15, 0]; % 1 triangle per row, position in m
tri_len = 0.01; % m
title_str = {'Loose Sand', 'Compacted Sand', 'Sieved Sand'};
title_str_start = [0, 146, 341]; 
scale_pos = [-0.1, 0.15, 0.1, 0.01]; % [x, y, width, height], in world coords
scale_label = '10 cm';
count_pos = [0.3, 0.15];

% define internal parameters
max_mov_dim = [1920, 1080];

tri_color = 'red';
tri_opacity = 0.7;

title_size = 72;
title_color = 'red';
title_box_color = 'white';
title_box_opacity = 1;

scale_color = 'red';
scale_opacity = 1;
scale_text_size = 72;
scale_text_color = 'red';
scale_box_color = 'white';
scale_box_opacity = 1;

count_size = 72;
count_color = 'red';
count_box_color = 'white';
count_box_opac = 1;

% parse inputs
if nargin < 3;
    show_frame = [];
end
if isempty(show_frame)
    make_movie = true;
else
    make_movie = false;
end

% get netcdf ids
ncid = netcdf.open(input_file, 'NOWRITE');
intens_id = netcdf.inqVarID(ncid, 'intensity');

% get coordinate vectors
x = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'x'));
y = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'y'));
step = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'step'));

% prepare movie object: 16-bit grayscale Motion JPEG 2000, lossless compression
if make_movie
    movie_writer = VideoWriter(output_file, 'Archival');
    movie_writer.MJ2BitDepth = 16;
    movie_writer.FrameRate = 10; % frames/second
    movie_writer.open();
end

% loop: read data, annotate images, create frames
for i = 1:numel(step)

    if ~make_movie && step(i) ~= show_frame
        continue
    end

    % update user
    fprintf('step: %i\n', step(i));
 
    % read intensity image
    intens = netcdf.getVar(ncid, intens_id, [0, 0, i-1], [numel(x), numel(y), 1])';
    
    % resize and reshape image as needed
    [frame, xf, yf] = yalebox_movie_aux_resize(intens, x, y, max_mov_dim);
    [frame, yf] = yalebox_movie_aux_flip(frame, yf);
    
    % add annotations (triangles, scale, title, counter)
    frame = yalebox_movie_aux_spoint(frame, xf, yf, tri_tip, tri_len, ...
        tri_color, tri_opacity);
    
    title_ind = find(step(i) >= title_str_start , 1, 'last');    
    frame = yalebox_movie_aux_title(frame, title_str{title_ind}, title_size, ...
        title_color, title_box_color, title_box_opacity);
    
    frame = yalebox_movie_aux_scalebar(frame, xf, yf, scale_label, scale_pos, ...
        scale_color, scale_opacity, scale_text_size, scale_text_color, ...
        scale_box_color, scale_box_opacity);
    
    frame = yalebox_movie_aux_counter(frame, xf, yf, step(i), count_pos, ...
        count_size, count_color, count_box_color, count_box_opac);

    % add frame to movie
    if make_movie
        writeVideo(movie_writer, im2frame(double(frame)));
    end
    
end

% finalize
netcdf.close(ncid);
if make_movie
    movie_writer.close();
else
    imshow(frame);
end