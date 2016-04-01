function [] = piv_check_manual(image_file, displ_file, step, result_file)
%
% Estimate displacements manually by selecting matching points in an image pair.
% Writes results to a csv file with columns: row_ti, column_ti, row_tf,
% column_tf, row_tm, column_tm, u_chk_tm, v_chk_tm, u_piv_tm, v_piv_tm.
%
% Arguments:
% 
% image_file = String, path to netCDF file containing PIV input images, as
%   produced by prep_series.m
%
% displ_file = String, path to netCDF file containing PIV output data, as
%   produced by piv_series.m
%
% step = Scalar, step to be checked, this should be a valid time step in
%   displ_file, and so will probably at the midpoint between two images (e.g.
%   11.5)
%
% result_file = String, path to the output CSV file, the program will prompt
%   before overwriting
% %

%% init

% define constants
default_num_pts = 20;

% sanity check
validateattributes(image_file,  {'char'}, {'vector'}, mfilename, 'image_file');
validateattributes(displ_file,  {'char'}, {'vector'}, mfilename, 'displ_file');
validateattributes(result_file, {'char'}, {'vector'}, mfilename, 'result_file');
validateattributes(step, {'numeric'}, {'scalar'}, mfilename, 'step');




%% extract the input data

% find the index of the displacements and images
displ_step = ncread(displ_file, 'step');
image_step = ncread(image_file, 'step');

displ_index_tm = find(displ_step == step);
assert(numel(displ_index_tm) == 1);

image_index_ti = displ_index_tm;
image_index_tf = displ_index_tm+1;
assert(image_step(image_index_ti) == floor(step));
assert(image_step(image_index_tf) == ceil(step));

% extract image, displacement, and coordinate data, transposing as needed
image_ini = double(ncread(image_file, 'intensity', [1, 1, image_index_ti], [inf, inf, 1]))';
image_fin = double(ncread(image_file, 'intensity', [1, 1, image_index_tf], [inf, inf, 1]))';
image_x = double(ncread(image_file, 'x'));
image_y = double(ncread(image_file, 'y'));

displ_u = double(ncread(displ_file, 'u', [1, 1, displ_index_tm], [inf, inf, 1]))';
displ_v = double(ncread(displ_file, 'v', [1, 1, displ_index_tm], [inf, inf, 1]))';
displ_mag = sqrt(displ_u.^2+displ_v.^2);
displ_x = double(ncread(displ_file, 'x'));
displ_y = double(ncread(displ_file, 'y'));

%% create the GUI

% create full-screen figure
f = figure;
f.Units = 'Normalized';
f.Position = [0, 0, 1, 1];

% create struct for shared data
f.UserData = struct(...
    'num_pts', []);

% plot velocity magnitude
ax_displ = axes('Position', [0.05, 0.7, 0.8, 0.25]);
imagesc(displ_x, displ_y, displ_mag, 'AlphaData', ~isnan(displ_mag));
ax_displ.YDir = 'Normal';
ax_displ.Color = 0.9*[1 1 1];
title('Displacement Magnitude')
axis equal tight
colorbar

% create axes for ini and fin
ax_ini = axes('Position', [0.05, 0.05, 0.35, 0.5]);
title('Initial Image');

ax_fin = axes('Position', [0.45, 0.05, 0.35, 0.5]);
title('Final Image');

% create num_pts control
text_num_pts = uicontrol('Style', 'text');
text_num_pts.Units = 'Normalized';
text_num_pts.Position = [0.85, 0.75, 0.1, 0.05];
text_num_pts.BackgroundColor = [1 1 1];
text_num_pts.String = 'Number of random test points';

edit_num_pts = uicontrol('Style', 'edit');
edit_num_pts.Units = 'Normalized';
edit_num_pts.Position = [0.87, 0.69, 0.05, 0.05];
edit_num_pts.BackgroundColor = [1 1 1];
edit_num_pts.String = num2str(default_num_pts);
edit_num_pts.Tag = 'edit_num_pts';
edit_num_pts.Callback = @set_num_pts;
edit_num_pts.Enable = 'on';

% NOTE: add button to update points

button_start_analysis = uicontrol('Style', 'pushbutton');
button_start_analysis.Units = 'Normalized';
button_start_analysis.Position = [0.85, 0.6, 0.1, 0.05];
button_start_analysis.String = 'Start Analysis';
button_start_analysis.Tag = 'button_start_analysis';
button_start_analysis.Callback = @start_analysis;
button_start_analysis.Enable = 'on';

% UserData contains the index of the current point
text_analysis_pt = uicontrol('Style', 'text');
text_analysis_pt.Units = 'Normalized';
text_analysis_pt.Position = [0.85, 0.5, 0.1, 0.05];
text_analysis_pt.BackgroundColor = [1 1 1];
text_analysis_pt.String = sprintf('Point 0 of %s', edit_num_pts.String);
text_analysis_pt.Tag = 'text_analysis_pt';
text_analysis_pt.UserData = 0;

% UserData contains the location of each control point as a num_pts*2 array:
% ...with a row [x, y] for each point
button_next_pt = uicontrol('Style', 'pushbutton');
button_next_pt.Units = 'Normalized';
button_next_pt.Position = [0.85, 0.45, 0.1, 0.05];
button_next_pt.String = 'Analyze Next Point';
button_next_pt.Tag = 'button_next_pt';
button_next_pt.Callback = {@next_pt, ax_ini, ax_fin, image_x, image_y, ...
                            image_ini, image_fin};
button_next_pt.Enable = 'off';

%% Initialize the GUI


end

function set_num_pts(hui, ~)
% Set the number of control points 

num_pts = str2double(hui.String);

if isempty(num_pts) || isnan(num_pts) || num_pts <= 0 || round(num_pts) ~= num_pts
    warndlg('Number of points must be a positive integer')
    return
else
    hui.Parent.UserData.num_pts = num_pts;
end

end



