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
edit_clim_min.Enable = 'on';

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
edit_clim_max.Enable = 'on';

text_num_pts = uicontrol('Style', 'text');
text_num_pts.Units = 'Normalized';
text_num_pts.Position = [0.85, 0.75, 0.1, 0.05];
text_num_pts.BackgroundColor = [1 1 1];
text_num_pts.String = 'Number of random test points';

% UserData contains input data for all control points as a num_pts*4 array:
% ... with a row [x, y, u, v] for each point
edit_num_pts = uicontrol('Style', 'edit');
edit_num_pts.Units = 'Normalized';
edit_num_pts.Position = [0.87, 0.69, 0.05, 0.05];
edit_num_pts.BackgroundColor = [1 1 1];
edit_num_pts.String = '20';
edit_num_pts.Tag = 'edit_num_pts';
edit_num_pts.Callback = {@update_test_pts, ax_displ, displ_x, displ_y, ...
                         ~isnan(displ_u), displ_u, displ_v};
edit_num_pts.Enable = 'on';

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

function next_pt(hObject, ~, ax_ini, ax_fin, xx, yy, ini, fin)
%
% Perform manual check for next point
%
% Arguments:
%
%   ~, ~ = unused, MATLAB GUI required arguments
%
%   xx, yy = Vector, coordinate vectors for ini and fin
%
%   ax_ini, ax_fin = Axes objects for initial and final image plots
%
%   ini, fin = 2D matrix, initial and final images
% %

% NOTE: estimated points are troublingly bad

% define parameters
image_x_dim = 0.02; % dimensions of image plots in world coordinates
image_y_dim = 0.02;

% get point index, update GUI
h = findobj('Tag', 'text_analysis_pt');
prev_pt_index = h.UserData;
pt_index = prev_pt_index+1;
h.UserData = pt_index;
h.String = num2str(pt_index);

% get control point location and PIV displacement
h = findobj('Tag', 'edit_num_pts');
pts = h.UserData;
num_pts = size(pts, 1);

if num_pts == 0;
    h_gui = gcf;
    errordlg('No test points are available - did you forget to generate them?');
    close(h_gui);
    return
end

x_pt = pts(pt_index, 1);
y_pt = pts(pt_index, 2);
u_pt_piv = pts(pt_index, 3);
v_pt_piv = pts(pt_index, 4);

% initialize output data the first time
if isempty(hObject.UserData)
    hObject.UserData = nan(num_pts, 2);
end

% get the estimated location of the control point at initial and final times
if isnan(hObject.UserData(pt_index, 1))
    x_pt_ini = x_pt-0.5*u_pt_piv;
    y_pt_ini = y_pt-0.5*v_pt_piv;
    x_pt_fin = x_pt+0.5*u_pt_piv;
    y_pt_fin = y_pt+0.5*v_pt_piv;
    
else
    x_pt_ini = x_pt-0.5*hObject.UserData(pt_index, 1);
    y_pt_ini = y_pt-0.5*hObject.UserData(pt_index, 2);
    x_pt_fin = x_pt+0.5*hObject.UserData(pt_index, 1);
    y_pt_fin = y_pt+0.5*hObject.UserData(pt_index, 2);

end
    
% plot estimated location of control points
axes(ax_ini);
imagesc(xx, yy, ini);
hold on
pt_ini = impoint(gca, x_pt_ini, y_pt_ini);
pt_ini.setColor('k');
xlim(x_pt_ini+image_x_dim*[-0.5, 0.5]);
ylim(y_pt_ini+image_y_dim*[-0.5, 0.5]);
hold off

axes(ax_fin);
imagesc(xx, yy, fin);
hold on
pt_fin = impoint(gca, x_pt_fin, y_pt_fin);
pt_fin.setColor('k');
xlim(x_pt_fin+image_x_dim*[-0.5, 0.5]);
ylim(y_pt_fin+image_y_dim*[-0.5, 0.5]);
hold off

% record control point location in UserData
hObject.UserData(pt_index, :) = pt_fin.getPosition()-pt_ini.getPosition();

disp(hObject.UserData);


function start_analysis(~, ~)
% function start_analysis(~, ~)
%
% Disables controls for initialization, enables controls for analysis
%
% Arguments:
%   ~ = unused, MATLAB GUI required arguments
% %

% disable initialization controls
h = findobj('Tag', 'edit_clim_min');         h.Enable = 'off';
h = findobj('Tag', 'edit_clim_max');         h.Enable = 'off';
h = findobj('Tag', 'edit_num_pts');          h.Enable = 'off';
h = findobj('Tag', 'button_start_analysis'); h.Enable = 'off';

% enable analysis controls
h = findobj('Tag', 'button_next_pt'); h.Enable = 'on';

function update_test_pts(hObject, ~, ax, xx, yy, roi, uu, vv)
%
% Generates a new batch of random test points, plots them on the displ axis, and
% returns their location in the .UserData handle property. All new points must
% lie with the ROI.
% 
% Arguments:
%
%   hObject = uicontrol object handle
%
%   ~ = unused, MATLAB GUI required arguments
%
%   ax = Axes handle for plot to be modified
%
%   xx, yy = Vector, coordinate vectors
%
%   roi = 2D, logical flags indicating points in the ROI (1) and not (0)
%
%   uu, vv = 2D, estimated displacement components at midpoint time
%
% Results are save in hObject.UserData as [x_pt, y_pt, u_pt, v_pt]
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

% create interpolant to estimate displacement at test points for midpoint time
[ynorm, xnorm] = ndgrid(linspace(0, 1, size(uu,1)), linspace(0, 1, size(uu,2)));
uu = flipud(uu);
vv = flipud(vv);
roi = logical(roi);
interpolant = scatteredInterpolant(ynorm(roi), xnorm(roi), uu(roi), 'linear');
u_pts = interpolant(x_pts, y_pts);
interpolant.Values = vv(roi);
v_pts = interpolant(x_pts, y_pts);

% convert normalized coordinates to world coordinates
x_pts = x_pts*range(xx)+min(xx);
y_pts = y_pts*range(yy)+min(yy);

% return results using UserData 
hObject.UserData = [x_pts, y_pts, u_pts, v_pts];

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


