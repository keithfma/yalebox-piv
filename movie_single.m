function [] = movie_single(prm, movie_type, show_frame)
% function [] = movie_single(prm, movie_type, show_frame)
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
% movie_type = String, select movie type, valid options are: color, streak
%
% show_frame = OPTIONAL, show a single frame at step == show_frame, do not
%   process other frames or make a movie, used for testing parameter
%   values.
% %

% parse inputs
make_movie = nargin < 3 || isempty(show_frame);
assert(ismember(movie_type, {'color', 'streak'}));
switch movie_type
    case 'color'        
        read_image = @read_image_color;
        output_file_matlab = [prm.output_stub '_color.mj2'];
        output_file_ffmpeg = [prm.output_stub, '_color.mp4'];
        output_file_ffmpeg_small = [prm.output_stub, '_color_small.mp4'];
        
    case 'streak'
        read_image = @read_image_gray;
        output_file_matlab = [prm.output_stub '_streak.mj2'];
        output_file_ffmpeg = [prm.output_stub, '_streak.mp4'];
        output_file_ffmpeg_small = [prm.output_stub, 'streak_small.mp4'];
        img_prev = 0;
        
end

% open netCDF file and get coordinate vars
ncid = netcdf.open(prm.input_file, 'NOWRITE');
x = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'x'));
y = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'y'));
step = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'step'));

% prepare movie object...format = motion JPEG 2000, lossless compression
if make_movie
    movie_writer = VideoWriter(output_file_matlab, 'Archival');
    movie_writer.FrameRate = prm.frame_rate; % frames/second
    movie_writer.open();
end

% loop: read data, annotate images, create frames
for i = 1:numel(step)
    
    if ~make_movie && step(i) ~= show_frame
        continue
    end
    
    % update user
    fprintf('step: %i\n', step(i));
    
    % read image
    img = read_image(ncid, numel(x), numel(y), i);
    
    % resize and reshape image as needed
    [img, xf, yf] = movie_frame_resize(img, x, y, prm.max_dim);
    [img, yf] = movie_frame_flip(img, yf);
    
    % streak only: apply threshold and decay
    if strcmp(movie_type, 'streak')        
        img(img<prm.streak_threshold) = 0;
        img(img>=prm.streak_threshold) = 1;
        img = img + img_prev*prm.streak_memory;
        img_prev = img; 
        img(img > 1) = 1;
    end
    
    % add annotations (triangles, scale, title, counter)
    frame = movie_frame_spoint(img, xf, yf, prm.tri_tip, prm.tri_len, ...
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
        writeVideo(movie_writer, im2frame(frame));
    end
    
end

% finalize
netcdf.close(ncid);

if ~make_movie
    % display processed frame
    figure
    imshow(frame);
else
    % finish movie and convert movie to better format using ffmpeg
    % ...MATLAB for Linux video support is limited, so I use ffmpeg
    % ...ffmpeg commands lines modified from https://trac.ffmpeg.org/wiki/Encode/H.264
    movie_writer.close();
    

    
    cmd_ffmpeg = sprintf('ffmpeg -i %s -c:v libx264 -preset veryslow -crf 18  -aspect:v %f -pix_fmt yuv420p %s', ...
        output_file_matlab, size(frame,2)/size(frame,1), output_file_ffmpeg);
    cmd_ffmpeg_small = sprintf('ffmpeg -i %s -c:v libx264 -preset veryslow -crf 28  -aspect:v %f -pix_fmt yuv420p %s', ...
        output_file_matlab, size(frame,2)/size(frame,1), output_file_ffmpeg_small);
    
    system(cmd_ffmpeg);
    system(cmd_ffmpeg_small);
    
    fprintf('%s: Complete\n', mfilename);
    fprintf('Generated the following files:\n');
    fprintf('-original (MATLAB): %s\n', output_file_matlab);
    fprintf('-high quality(FFMPEG): %s\n', output_file_ffmpeg);
    fprintf('-small (FFMPEG): %s\n', output_file_ffmpeg_small);
end


function im = read_image_gray(ncid, nx, ny, step)
% Read grayscale image from netCDF variable 'img'

im = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'img'), [0, 0, step-1], [ny, nx, 1]);
im = double(im);


function im = read_image_color(ncid, nx, ny, step)
% Read color image from netCDF variable 'img_rgb'

im = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'img_rgb'), [0, 0, 0, step-1], [ny, nx, 3, 1]);
mask_auto = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'mask_auto'), [0, 0, step-1], [ny, nx, 1]);
mask_manual = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'mask_manual'));
mask = mask_auto & mask_manual;
im(repmat(~mask, [1, 1, 3])) = 0;