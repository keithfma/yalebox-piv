function [] = movie_single(prm, show_frame)
% function [] = movie_single(prm, show_frame)
% 
% Generate video for a single view of a yalebox experiment from the output of
% the pre-processing routine prep_series(). The (many) parameters below control
% the content, size, annotation, memory, etc for the video. All processing is
% done in MATLAB, which restricts the type of the output movie file. This can be
% converted easily after the fact with ffmpeg, and the program prints the
% command line needed to do so.
%
% ----- Source parameters -----
% 
% prm.input_file = String, path to input netCDF file as produced by
%   prep_series()
%
% prm.output_stub = String, path to output video file, without file extension
%
% prm.input_var = String, name of the netCDF variable that contains the input
%   images, valid selections are: img, img_raw
% 
% ----- Video parameters -----
%
% prm.frame_rate = Scalar integer, frames/sec in output video
%
% prm.max_dim = Vector, length == 2, integer, [row, col] maximum dimensions of
%   the output video, image is resized to fit in this bounding box
%
% prm.threshold = Scalar, threshold value applied to image, pixels <
%   this value are set to zero
%
% prm.memory = Scalar in the range [0, 1], controls the decay of the previous
%   image, where 1 adds the previous to the current image, and 0 adds nothing from
%   the previous image
%
% ----- S-point triangle annotation parameters -----
%
% prm.tri_tip = Vector, size = [number of triangles, 2], tip location of s-point
%   triangles in [x, y] world coordinates, with 1 triangle in each row 
%
% prm.tri_len = Scalar, side length for all (equilateral) triangles in world
%   coordinates 
%
% prm.tri_color = String, color definition MATLAB can understand, e.g. 'red'
%
% prm.tri_opacity = Scalar range [0, 1], "alpha" for the triangle patch
%
% ----- Title annotation parameters -----
%
% prm.title_str = Cell array of strings, title string for each experiment
%   segment
%
% prm.title_str_start = Vector, starting index for each experiment segment 
%
% prm.title_size = Scalar, integer, title font size in points
%
% prm.title_color = String, color definition MATLAB can understand, e.g. 'red'
%
% prm.title_box_color = String, color definition MATLAB can understand
%
% prm.title_box_opacity = Scalar, range [0, 1], "alpha" for the title box
% 
% ----- Scalebar annotation parameters ----- 
%
% prm.scale_pos = Vector, scalebar location and size as a MATLAB position vector
%   in world units, [left, bottom, width, height]
%
% prm.scale_label = String, scalebar label text
%
% prm.scale_color = String, color definition MATLAB can understand, e.g. 'red'
%
% prm.scale_opacity = Scalar, range [0, 1], "alpha" for the scalebar
%
% prm.scale_text_size = Scalar, integer, scale label font size in points
%
% prm.scale_text_color = String, color definition MATLAB can understand
%
% prm.scale_box_color = String, color definition MATLAB can understand
%
% prm.scale_box_opacity = Scalar, range [0, 1], "alpha" for the scale label box
%
% ----- Counter annotation parameters ----- 
%
% prm.count_pos = Vector, counter position in [x, y] world units
%
% prm.count_size = Scalar, integer, counter font size in points
%
% prm.count_color = String, color definition MATLAB can understand, e.g. 'red'
%
% prm.count_box_color = String, color definition MATLAB can understand, e.g. 'red'
%
% prm.count_box_opac = Scalar, range [0, 1], "alpha" for the counter text box
% 
% ----- Other -----
%
% show_frame = OPTIONAL, show a single frame at step == show_frame, do not
%   process other frames or make a movie, used for testing parameter
%   values.
% %

% parse inputs
make_movie = nargin < 4 || isempty(show_frame);

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
    movie_writer.FrameRate = prm.frame_rate; % frames/second
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
    intens(intens<prm.threshold) = 0;
    intens = intens + intens_prev*prm.memory;
    intens_prev = intens; % prep for next frame
    
    % resize and reshape image as needed
    [frame, xf, yf] = movie_frame_resize(intens, x, y, prm.max_dim);
    [frame, yf] = movie_frame_flip(frame, yf);
    
    % add annotations (triangles, scale, title, counter)
    frame = movie_frame_spoint(frame, xf, yf, prm.tri_tip, prm.tri_len, ...
        prm.tri_color, prm.tri_opacity);
    
    title_ind = find(step(i) >= prm.title_str_start , 1, 'last');    
    frame = movie_frame_title(frame, prm.title_str{title_ind}, ...
        prm.title_size, prm.title_color, prm.title_box_color, prm.title_box_opacity);
    
    frame = movie_frame_scalebar(frame, xf, yf, prm.scale_label, ...
        prm.scale_pos, prm.scale_color, prm.scale_opacity, prm.scale_text_size, ...
        prm.scale_text_color, prm.scale_box_color, prm.scale_box_opacity);
    
    frame = movie_frame_counter(frame, xf, yf, step(i), prm.count_pos, ...
        prm.count_size, prm.count_color, prm.count_box_color, prm.count_box_opac);

    % add frame to movie
    if make_movie
        writeVideo(movie_writer, im2frame(double(frame)));
    end
    
end

% finalize
netcdf.close(ncid);
if make_movie
    movie_writer.close();
    fprintf('Run the following command to re-encode as an MP4 video:\n');
    fprintf('\tavconv <original video> -c:v libx264 -aspect:v %f <new mp4 video>\n', ...
        size(frame, 2)/size(frame, 1));    
else
    imshow(frame);
end

