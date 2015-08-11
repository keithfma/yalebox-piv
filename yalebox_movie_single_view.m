function [] = yalebox_movie_single_view(input_file, output_file)
%
% input_file = String, path to input netCDF file as produced by
%   yalebox_prep_create_input().
%
% output_file = String path to output video file, without file extension
%
% Keith Ma, August 2015

% debug arguments
tri_tip_xy = [0, 0; -0.15, 0]; % 1 triangle per row, position in m
tri_side_len = 0.01; % m

title_str = 'Party Time';

% internal parameters
tri_color = 'red';
tri_opacity = 0.7;
title_font_size = 20;
title_font_color = 'red';

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

% prepare annotation arguments
tri_pos = get_tri_pos(tri_tip_xy, tri_side_len, x, y);



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
        
    % add s-point triangles
    frame = insertShape(intens, 'FilledPolygon', tri_pos, ...
        'Color', tri_color, ...
        'Opacity', tri_opacity);
    
    frame = insertText(frame, [nx/2, 1], title_str, ...
        'FontSize', title_font_size, ...
        'TextColor', title_font_color, ...
        'BoxOpacity', 0, ...
        'AnchorPoint', 'CenterTop');
        
    % create frame    
end

% finalize
netcdf.close(ncid);
% movie_writer.close();

keyboard

% helper functions --------------------------------------------------------

function pos = get_tri_pos(tip, len, x, y)
%
% Return position argument used by insertShape to draw s-point triangles
%

% init
nx = numel(x);
ny = numel(y);
nt = size(tip, 1); % 1 triangle per row
dx = abs(x(1)-x(2));
dy = abs(y(1)-y(2));
pos = nan(nt, 6); 

for i = 1:nt
   rr = interp1(y, 1:ny, tip(i,2), 'linear', 'extrap'); 
   cc = interp1(x, 1:nx, tip(i,1), 'linear', 'extrap');
      
   if abs(rr-ny) < abs(rr) % tip at top
       pos(i,:) = [cc,          rr,                  ...
                   cc+len/2/dx, rr-sqrt(3)/2*len/dy, ...
                   cc-len/2/dx, rr-sqrt(3)/2*len/dy];       
   
   else % tip at bottom
       pos(i,:) = [cc,          rr,                  ...
                   cc+len/2/dx, rr+sqrt(3)/2*len/dy, ...
                   cc-len/2/dx, rr+sqrt(3)/2*len/dy];
   end   
end