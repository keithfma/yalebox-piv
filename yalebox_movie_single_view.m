function [] = yalebox_movie_single_view(input_file, output_file)
%
% input_file = String, path to input netCDF file as produced by
%   yalebox_prep_create_input().
%
% output_file = String path to output video file, without file extension
%
% Keith Ma, August 2015

% define parameters
len_tri = 0.01; % m
color_tri = 'yellow';
opacity = 0.7;


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

% get shape annotations

spt_r = interp1(y, 1:ny, 0, 'linear', 'extrap');
spt_c = interp1(x, 1:nx, 0, 'linear', 'extrap');
if abs(spt_r-ny) < abs(spt_r) % spt at image top
    pos_tri = [spt_c, spt_r, ...
               spt_c+len_tri/2/dx, spt_r-sqrt(3)/2*len_tri/dy, ...
               spt_c-len_tri/2/dx, spt_r-sqrt(3)/2*len_tri/dy];
else % spt at image bottom
end

shape_pos = {pos_tri};
shape_color = {color_tri};



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
    intens = insertShape(intens, 'FilledPolygon', shape_pos, ...
        'Color', shape_color, 'Opacity', opacity);
    
    
    % create frame    
end

% finalize
netcdf.close(ncid);
% movie_writer.close();

keyboard