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

% control for velocity color scale (two edit boxes and a reset button)
clim = caxis(ax_displ);

text_clim_min = uicontrol('Style', 'text');
text_clim_min.Units = 'Normalized';
text_clim_min.Position = [0.85, 0.9, 0.05, 0.05];
text_clim_min.BackgroundColor = [1 1 1];
text_clim_min.String = 'Min [mm/step]';

edit_clim_min = uicontrol('Style', 'edit');
edit_clim_min.Units = 'Normalized';
edit_clim_min.Position = [0.85, 0.85, 0.05, 0.05];
edit_clim_min.BackgroundColor = [1 1 1];
edit_clim_min.String = num2str(clim(1)*1000);
edit_clim_min.Tag = 'edit_clim_min';
edit_clim_min.Callback = {@update_displ_clim, ax_displ};

text_clim_max = uicontrol('Style', 'text');
text_clim_max.Units = 'Normalized';
text_clim_max.Position = [0.92, 0.9, 0.05, 0.05];
text_clim_max.BackgroundColor = [1 1 1];
text_clim_max.String = 'Max [mm/step]';

edit_clim_max = uicontrol('Style', 'edit');
edit_clim_max.Units = 'Normalized';
edit_clim_max.Position = [0.92, 0.85, 0.05, 0.05];
edit_clim_max.BackgroundColor = [1 1 1];
edit_clim_max.String = num2str(clim(2)*1000);
edit_clim_max.Tag = 'edit_clim_max';
edit_clim_max.Callback = {@update_displ_clim, ax_displ};

text_num_pts = uicontrol('Style', 'text');
text_num_pts.Units = 'Normalized';
text_num_pts.Position = [0.85, 0.75, 0.1, 0.05];
text_num_pts.BackgroundColor = [1 1 1];
text_num_pts.String = 'Number of random test points';

edit_num_pts = uicontrol('Style', 'edit');
edit_num_pts.Units = 'Normalized';
edit_num_pts.Position = [0.87, 0.69, 0.05, 0.05];
edit_num_pts.BackgroundColor = [1 1 1];
edit_num_pts.String = '20';
edit_num_pts.Tag = 'edit_num_pts';
edit_num_pts.Callback = {@update_test_pts, displ_x, displ_y, ~isnan(displ_u), ax_displ};

% GUI Functions ----------------------------------------------------------------

function update_test_pts(hObject, ~, roi_x, roi_y, roi, ax)
%
% Generates a new batch of random test points, plots them on the displ axis, and
% returns their location in the .UserData handle property. All new points must
% lie with the ROI.
% 
% Arguments:
%   hObject = uicontrol object handle
%   ~ = unused, MATLAB GUI required arguments
%   x_roi, y_roi = Vector, coordinate vectors for the ROI
%   roi = 2D, logical flags indicating points in the ROI (1) and not (0)
%   ax = Axes handle for plot to be modified
% % 

% get number of test points
h_num_pts = findobj('Tag', 'edit_num_pts');
num_pts = str2double(h_num_pts.String);

if isempty(num_pts) || isnan(num_pts) || num_pts <= 0
    warndlg('invalid value for number of test points');
    return
end

% create interpolant to check if test points lie in the ROI, use increasing normalized coordinates
[ynorm, xnorm] = ndgrid(linspace(0, 1, size(roi,1)), linspace(0, 1, size(roi,2)));
roi = double(flipud(roi));
interpolant = griddedInterpolant(ynorm, xnorm, roi, 'nearest');

% generate test points, use normalized coordinates, keep only if in ROI
x_pts = nan(num_pts, 1);
y_pts = nan(num_pts, 1);
for ii = 1:num_pts
    while isnan(x_pts(ii))
       x_candidate = rand(1);
       y_candidate = rand(1);
       if interpolant(y_candidate, x_candidate) == 1
           x_pts(ii) = x_candidate;
           y_pts(ii) = y_candidate;
           break
       end
    end
end

% convert normalized coordinates to world coordinates
x_pts = x_pts*range(roi_x)+min(roi_x);
y_pts = y_pts*range(roi_y)+min(roi_y);

% return results using UserData 
hObject.UserData = [x_pts, y_pts];

% clear previous points from axes, then plot new points
axes(ax);
h_prev = findobj(gca, 'Type', 'Line');
delete(h_prev);

hold on
plot(x_pts, y_pts, 'xk');


function update_displ_clim(~, ~, ax)
% function update_displ_clim(~, ~, ax)
%
% Updates displacement plot color axis using the values in the two edit boxes.
%
% Arguments:
%   ~, ~ = unused, MATLAB GUI required arguments
%   ax = Axes object for the displacement plot
% %

% get requested caxis
hmin = findobj('Tag', 'edit_clim_min');
hmax = findobj('Tag', 'edit_clim_max');
clim_min = str2double(hmin.String);
clim_max = str2double(hmax.String);

% sanity check
if isempty(clim_min) || isnan(clim_min)
    warndlg('invalid value for minimum color in displacement plot');
    return
end
if isempty(clim_max) || isnan(clim_max)
    warndlg('invalid value for minimum color in displacement plot');
    return
end
if clim_min >= clim_max
    warndlg('minimum value is greater than maximum value in displacement plot');
    return
end

% update colors
caxis(ax, [clim_min, clim_max]/1000);







% % 
% % % % Edit: number of test points
% % % text_num_pts = uicontrol('Style', 'text');
% % % text_num_pts.Units = 'Normalized';
% % % text_num_pts.Position = [control_panel_left, 0.9, control_width, control_height];
% % % text_num_pts.String = 'Number of control points:';
% % % text_num_pts.BackgroundColor = [1 1 1];
% % % 
% % % edit_num_pts = uicontrol('Style', 'edit');
% % % edit_num_pts.Units = 'Normalized';
% % % edit_num_pts.Position = [control_panel_left, control_panel_top-1.1*control_height, control_width, control_height];
% % % edit_num_pts.String = '?';
% % % edit_num_pts.Tag = 'edit_num_pts';
% % % edit_num_pts.Callback = @generate_pts;
% % 
% % % % Button: recompute color limits
% % % button_reset_clim = uicontrol('Style', 'pushbutton');
% % % button_reset_clim.Units = 'Normalized';
% % % button_reset_clim.Position = [control_panel_left, 0.9, button_width, button_height];
% % % button_reset_clim.String = 'Reset Colors';
% % % button_reset_clim.Callback = {@callback_button_reset_clim, ax_mag, displ_x, displ_y, displ_mag};
% % % 
% % % % Button: add point
% % % button_add_point = uicontrol('Style', 'pushbutton');
% % % button_add_point.Units = 'Normalized';
% % % button_add_point.Position = [control_panel_left, 0.8, button_width, button_height];
% % % button_add_point.String = 'Add Point';
% % % button_add_point.Callback = {@callback_button_add_point, ax_mag, ax_ini, ax_fin, displ_x, displ_y, displ_u, displ_v};
% % % button_add_point.UserData = 0;
% % 
% % % Button: delete point
% % 
% % end
% % 
% % % GUI Functions ----------------------------------------------------------------
% % 
% % function generate_pts(~, ~)
% % % function edit_num_pts_callback(hObject, ~)
% % %
% % % Updates the number of random points, then regenerates the points.
% % %
% % % Arguments:
% % %   hObject = Handle to self
% % %   ~ = unused, MATLAB GUI required arguments
% % % %
% % 
% % % get number of points 
% % h = findobj('Tag', 'edit_num_pts');
% % num_pts = str2double(h.String);
% % disp(num_pts)
% % if isempty(num_pts) || round(num_pts) ~= num_pts
% %     warndlg(sprintf('Expected an integer number of points, received "%s".', h.String));
% %     return
% % end
% % 
% % % generate random points
% % 
% % end
% % 
% % % function callback_button_reset_clim(~, ~, ax, cc, rr, zz)
% % % %
% % % % Resets the colors of axes "ax" to the [min, max] of the displayed area
% % % %
% % % % Arguments:
% % % %   ~, ~ = unused, MATLAB GUI required arguments
% % % %   ax = Axes object for plot to be rescaled
% % % %   cc, rr = Coordinate vectors (columns and rows) for plot to be rescaled
% % % %   zz = Gridded data for plot to be rescaled
% % % % %
% % % 
% % % visible_col_min = find(cc > ax.XLim(1), 1, 'first'); 
% % % visible_col_max = find(cc < ax.XLim(2), 1, 'last'); 
% % % visible_row_min = find(rr > ax.YLim(1), 1, 'first'); 
% % % visible_row_max = find(rr < ax.YLim(2), 1, 'last'); 
% % % visible_zz = zz(visible_row_min:visible_row_max, visible_col_min:visible_col_max);
% % % ax.CLim = [min(visible_zz(:)), max(visible_zz(:))];
% % % 
% % % end
% % % 
% % % function callback_button_add_point(hObject, ~, ax_mag, ax_ini, ax_fin, displ_x, displ_y, displ_u, displ_v, label)
% % % %
% % % % Add a new point to the analysis. This callback waits for the user to click a
% % % % point on the displacement magnitude axes, then creates interactive, editable,
% % % % impoint objects in the initial and final image axes. Clicking outside the
% % % % displacement magnitude axes generates a warning dialog box and does not add
% % % % any points.
% % % %
% % % % Arguments:
% % % %   hObject = Handle to self
% % % %   ~ = unused, MATLAB GUI required arguments
% % % %   ax_mag, ax_ini, ax_fin = Axes objects for displacement magnitude and initial
% % % %       and final images.
% % % %   displ_x, displ_y = Coordinate vectors for displacement grids
% % % %   displ_u, displ_v = Displacement grids
% % % % %
% % % 
% % % % select point from displacement magnitude axes
% % % axes(ax_mag)
% % % [x_tm, y_tm] = ginput(1);
% % % 
% % % % exit with warning if user selected a point outside the visible area
% % % if x_tm < ax_mag.XLim(1) || x_tm > ax_mag.XLim(2) || y_tm < ax_mag.YLim(1) || y_tm > ax_mag.YLim(2)
% % %     warndlg('Selected point is outside the displacement magnitude axes, skipping');
% % %     return
% % % end
% % % 
% % % % interpolate displacements at the selected point
% % % [displ_x_grid, displ_y_grid] = meshgrid(displ_x, displ_y);
% % % data = ~isnan(displ_u) & ~isnan(displ_v);
% % % interpolant = scatteredInterpolant(displ_x_grid(data), displ_y_grid(data), displ_u(data), 'linear');
% % % u_tm = interpolant(x_tm, y_tm);
% % % interpolant.Values = displ_v(data);
% % % v_tm = interpolant(x_tm, y_tm);
% % % 
% % % % get label and update counter
% % % label_str = num2str(hObject.UserData);
% % % 
% % % % plot point as label
% % % text(x_tm, y_tm, label_str);
% % % hObject.UserData = hObject.UserData+1; 
% % % 
% % % % create impoint objects 
% % % h_pt_ini = impoint(ax_ini, x_tm-0.5*u_tm, y_tm-0.5*v_tm);
% % % h_pt_ini.setColor('k');
% % % h_pt_ini.setString(label_str);
% % % 
% % % h_pt_fin = impoint(ax_fin, x_tm+0.5*u_tm, y_tm+0.5*v_tm);
% % % h_pt_fin.setColor('k');
% % % h_pt_fin.setString(label_str);
% % % 
% % % end
% % 
% % %% Compute and write results
% % 
% % %% Restart from results file
