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
%
% Plan:
% 
% GUI with velocity magnitude, ini and fin in subplots, linked axes
% Button: add point, use impoint() to select a point in both images
% Button: Compute results: gather data from all the points, compute displacements and write the output csv
% Button: Restart from results file: read a results file and add impoint objects to restart
% %

%% check for sane inputs

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
image_ini = ncread(image_file, 'intensity', [1, 1, image_index_ti], [inf, inf, 1])';
image_fin = ncread(image_file, 'intensity', [1, 1, image_index_tf], [inf, inf, 1])';
image_x = ncread(image_file, 'x');
image_y = ncread(image_file, 'y');

displ_u = ncread(displ_file, 'u', [1, 1, displ_index_tm], [inf, inf, 1])';
displ_v = ncread(displ_file, 'v', [1, 1, displ_index_tm], [inf, inf, 1])';
displ_mag = sqrt(displ_u.^2+displ_v.^2);
displ_x = ncread(displ_file, 'x');
displ_y = ncread(displ_file, 'y');

%% create the GUI

% parameters
plot_panel_width = 0.8;
control_panel_width = 1-plot_panel_width;
control_panel_left = plot_panel_width;
plot_vert_spc = 0.1;
plot_horiz_spc = 0.05;
ax_width = plot_panel_width-2*plot_horiz_spc;
ax_height = (1-4*plot_vert_spc)/3;
color_nodata = 0.9*[1 1 1];
button_width = 0.1;
button_height = 0.05; 

% create full-screen figure
f = figure;
f.Units = 'Normalized';
f.Position = [0, 0, 1, 1];

% plot velocity magnitude, ini, and fin
ax_fin = axes('Position', [plot_horiz_spc, plot_vert_spc              , ax_width, ax_height]);
ax_ini = axes('Position', [plot_horiz_spc, 2*plot_vert_spc+ax_height  , ax_width, ax_height]);
ax_mag = axes('Position', [plot_horiz_spc, 3*plot_vert_spc+2*ax_height, ax_width, ax_height]);

axes(ax_mag);
imagesc(displ_x, displ_y, displ_mag, 'AlphaData', ~isnan(displ_mag));
ax_mag.YDir = 'Normal';
ax_mag.Color = color_nodata;
axis equal tight
title('Displacement Magnitude')

axes(ax_ini)
imagesc(image_x, image_y, image_ini, 'AlphaData', image_ini ~= 0)
ax_ini.YDir = 'Normal';
ax_ini.Color = color_nodata;
axis equal tight
title('Initial Image');

axes(ax_fin)
imagesc(image_x, image_y, image_fin, 'AlphaData', image_fin ~= 0)
ax_fin.YDir = 'Normal';
ax_fin.Color = color_nodata;
axis equal tight
title('Final Image');

linkaxes([ax_mag, ax_ini, ax_fin], 'xy');

% Button: recompute color limits
button_reset_clim = uicontrol('Style', 'pushbutton');
button_reset_clim.Units = 'Normalized';
button_reset_clim.Position = [control_panel_left, 0.9, button_width, button_height];
button_reset_clim.String = 'Reset Colors';
button_reset_clim.Callback = {@callback_button_reset_clim, ax_mag, displ_x, displ_y, displ_mag};

keyboard

end

% GUI Functions ----------------------------------------------------------------

function callback_button_reset_clim(hObject, callbackdata, ax, cc, rr, zz)
%
% Resets the colors of axes "ax" to the [min, max] of the displayed area
%
% Arguments:
%   hObject, callbackdata = MATLAB GUI required arguments
%   ax = Axes object for plot to be rescaled
%   cc, rr = Coordinate vectors (columns and rows) for plot to be rescaled
%   zz = Gridded data for plot to be rescaled
% %

visible_col_min = find(cc > ax.XLim(1), 1, 'first'); 
visible_col_max = find(cc < ax.XLim(2), 1, 'last'); 
visible_row_min = find(rr > ax.YLim(1), 1, 'first'); 
visible_row_max = find(rr < ax.YLim(2), 1, 'last'); 
visible_zz = zz(visible_row_min:visible_row_max, visible_col_min:visible_col_max);
ax.CLim = [min(visible_zz(:)), max(visible_zz(:))];

end






%% Compute and write results

%% Restart from results file
