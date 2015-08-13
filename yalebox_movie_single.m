function [] = yalebox_movie_single(input_file, output_file, p, show_frame)
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

% parse inputs
if nargin < 4;
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
    movie_writer.FrameRate = p.frame_rate; % frames/second
    movie_writer.open();
end

% loop: read data, annotate images, create frames
intens_prev = zeros(numel(y), numel(x));
for i = 1:numel(step)

    if ~make_movie && step(i) ~= show_frame
        continue
    end

    % update user
    fprintf('step: %i\n', step(i));
 
    % read intensity image
    intens = netcdf.getVar(ncid, intens_id, [0, 0, i-1], [numel(x), numel(y), 1])';
    
    % apply threshold and decay
    intens(intens<p.threshold) = 0;
    intens = intens + intens_prev*p.decay_factor;
    intens_prev = intens; % prep for next frame
    
    % resize and reshape image as needed
    [frame, xf, yf] = yalebox_movie_aux_resize(intens, x, y, p.max_dim);
    [frame, yf] = yalebox_movie_aux_flip(frame, yf);
    
    % add annotations (triangles, scale, title, counter)
    frame = yalebox_movie_aux_spoint(frame, xf, yf, p.tri_tip, p.tri_len, ...
        p.tri_color, p.tri_opacity);
    
    title_ind = find(step(i) >= p.title_str_start , 1, 'last');    
    frame = yalebox_movie_aux_title(frame, p.title_str{title_ind}, ...
        p.title_size, p.title_color, p.title_box_color, p.title_box_opacity);
    
    frame = yalebox_movie_aux_scalebar(frame, xf, yf, p.scale_label, ...
        p.scale_pos, p.scale_color, p.scale_opacity, p.scale_text_size, ...
        p.scale_text_color, p.scale_box_color, p.scale_box_opacity);
    
    frame = yalebox_movie_aux_counter(frame, xf, yf, step(i), p.count_pos, ...
        p.count_size, p.count_color, p.count_box_color, p.count_box_opac);

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