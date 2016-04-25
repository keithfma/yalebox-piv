function [] = piv_check_manual(image_file, displ_file, step, result_prefix)
% function [] = piv_check_manual(image_file, displ_file, step, result_prefix)
%
% Estimate displacements manually by selecting matching points in an image pair.
% Writes results to a csv file with columns: row_ti, column_ti, row_tf,
% column_tf, row_tm, column_tm, u_chk_tm, v_chk_tm, u_piv_tm, v_piv_tm.
%
% Arguments:
% 
%   image_file = String, path to netCDF file containing PIV input images, as
%       produced by prep_series.m
%
%   displ_file = String, path to netCDF file containing PIV output data, as
%       produced by piv_series.m
%
%   step = Scalar, step to be checked, this should be a valid time step in
%       displ_file, and so will probably at the midpoint between two images
%       (e.g. 11.5)
%
%   result_prefix = String, prefix for results files, can include the path. The
%       program will prompt before overwriting
% %

%% init

% sanity check
validateattributes(image_file,  {'char'}, {'vector'}, mfilename, 'image_file');
validateattributes(displ_file,  {'char'}, {'vector'}, mfilename, 'displ_file');
validateattributes(result_prefix, {'char'}, {'vector'}, mfilename, 'result_prefix');
validateattributes(step, {'numeric'}, {'scalar'}, mfilename, 'step');
    
% create struct for shared data
share = struct(...
    'image_file', image_file, ...
    'displ_file', displ_file, ...
    'step', step, ...
    'result_nc', [result_prefix, '.nc'], ...
    'result_fig_loc', [result_prefix, '_loc.fig'], ...
    'result_fig_pdf', [result_prefix, '_pdf.fig'], ...
    'result_fig_vs', [result_prefix, '_vs.fig'], ...
    'ini', [], ...
    'fin', [], ...
    'image_x', [], ...
    'image_y', [], ...    
    'piv_x', [], ...
    'piv_y', [], ...
    'piv_u', [], ...
    'piv_v', [], ...
    'piv_m', [], ...
    'num_pts', 20, ...
    'ctrl_x', [], ...
    'ctrl_y', [], ...
    'ctrl_u', [], ...
    'ctrl_v', [], ...
    'ii', 0);

% protect against accidental overwriting
if exist(share.result_nc, 'file') == 2 ...
        || exist(share.result_fig_loc, 'file') == 2 ...
        || exist(share.result_fig_pdf, 'file') == 2 ...
        || exist(share.result_fig_vs, 'file') == 2
    
    but = questdlg('One or more results file exists', 'WARNING', 'OVERWRITE','EXIT','EXIT', 'EXIT');
    if strcmp(but, 'OVERWRITE') ~= 1
        return
    end
end

% find the index of the displacements and images
displ_step = ncread(displ_file, 'step');
displ_index_tm = find(displ_step == step);
assert(numel(displ_index_tm) == 1);

image_step = ncread(image_file, 'step');
image_index_ti = displ_index_tm;
image_index_tf = displ_index_tm+1;
assert(image_step(image_index_ti) == floor(step));
assert(image_step(image_index_tf) == ceil(step));

% extract image, displacement, and coordinate data, transposing as needed
share.ini = double(ncread(image_file, 'intensity', [1, 1, image_index_ti], [inf, inf, 1]))';
share.fin = double(ncread(image_file, 'intensity', [1, 1, image_index_tf], [inf, inf, 1]))';
share.image_x = double(ncread(image_file, 'x'));
share.image_y = double(ncread(image_file, 'y'));

share.piv_u = double(ncread(displ_file, 'u', [1, 1, displ_index_tm], [inf, inf, 1]))';
share.piv_v = double(ncread(displ_file, 'v', [1, 1, displ_index_tm], [inf, inf, 1]))';
share.piv_m = sqrt(share.piv_u.^2 + share.piv_v.^2);
share.piv_roi = double(~isnan(share.piv_m));
share.piv_x = double(ncread(displ_file, 'x'));
share.piv_y = double(ncread(displ_file, 'y'));

% create full-screen figure and attach shared data, invisible until init is completed
hf = figure;
hf.Units = 'Normalized';
hf.Position = [0, 0, 1, 1];
hf.UserData = share;
hf.Visible = 'off';

% create axes 
axes('Position', [0.05, 0.7, 0.8, 0.25]);
set(gca, 'Tag', 'ax_displ', 'NextPlot','add'); % NextPlot protects the Tag value
title('Displacement Magnitude')

axes('Position', [0.05, 0.05, 0.35, 0.5]);
imagesc(share.piv_x, share.piv_y, share.ini);
hpt = impoint(gca, [0, 0]);
hpt.setColor('k');
set(gca, 'Tag', 'ax_ini', 'NextPlot','add', 'YDir', 'Normal', 'UserData', hpt);
title('Initial Image');

axes('Position', [0.45, 0.05, 0.35, 0.5]);
imagesc(share.piv_x, share.piv_y, share.fin);
hpt = impoint(gca, [0, 0]);
hpt.setColor('k');
set(gca, 'Tag', 'ax_fin', 'NextPlot','add', 'YDir', 'Normal', 'UserData', hpt);
title('Final Image');

% create num_pts control
text_num_pts = uicontrol('Style', 'text');
text_num_pts.Units = 'Normalized';
text_num_pts.Position = [0.85, 0.85, 0.1, 0.05];
text_num_pts.BackgroundColor = [1 1 1];
text_num_pts.String = 'Number of random test points';

edit_num_pts = uicontrol('Style', 'edit');
edit_num_pts.Units = 'Normalized';
edit_num_pts.Position = [0.87, 0.79, 0.05, 0.05];
edit_num_pts.BackgroundColor = [1 1 1];
edit_num_pts.String = share.num_pts;
edit_num_pts.Tag = 'edit_num_pts';
edit_num_pts.Callback = @set_num_pts;
edit_num_pts.Enable = 'on';

but_get_pts = uicontrol('Style', 'pushbutton');
but_get_pts.Units = 'Normalized';
but_get_pts.Position = [0.85, 0.7, 0.1, 0.05];
but_get_pts.String = 'Get Control Points';
but_get_pts.Tag = 'but_get_pts';
but_get_pts.Callback = @get_ctrl_pts;
but_get_pts.Enable = 'on';

button_start_analysis = uicontrol('Style', 'pushbutton');
button_start_analysis.Units = 'Normalized';
button_start_analysis.Position = [0.85, 0.6, 0.1, 0.05];
button_start_analysis.String = 'Start Analysis';
button_start_analysis.Tag = 'button_start_analysis';
button_start_analysis.Callback = @start_analysis;
button_start_analysis.Enable = 'on';

text_curr_pt = uicontrol('Style', 'text');
text_curr_pt.Units = 'Normalized';
text_curr_pt.Position = [0.85, 0.5, 0.1, 0.05];
text_curr_pt.BackgroundColor = [1 1 1];
text_curr_pt.String = '';
text_curr_pt.Tag = 'text_curr_pt';

but_next_pt = uicontrol('Style', 'pushbutton');
but_next_pt.Units = 'Normalized';
but_next_pt.Position = [0.85, 0.45, 0.1, 0.05];
but_next_pt.String = 'Next Point';
but_next_pt.Tag = 'but_next_pt';
but_next_pt.Callback = {@next_pt, 1, 1};
but_next_pt.Enable = 'off';

but_prev_pt = uicontrol('Style', 'pushbutton');
but_prev_pt.Units = 'Normalized';
but_prev_pt.Position = [0.85, 0.35, 0.1, 0.05];
but_prev_pt.String = 'Previous Point';
but_prev_pt.Tag = 'but_prev_pt';
but_prev_pt.Callback = {@next_pt, 1, -1};
but_prev_pt.Enable = 'off';

but_delete = uicontrol('Style', 'pushbutton');
but_delete.Units = 'Normalized';
but_delete.Position = [0.85, 0.25, 0.1, 0.05];
but_delete.String = 'Delete Point';
but_delete.Tag = 'but_delete';
but_delete.Callback = @delete_pt;
but_delete.Enable = 'off';

but_done = uicontrol('Style', 'pushbutton');
but_done.Units = 'Normalized';
but_done.Position = [0.85, 0.15, 0.1, 0.05];
but_done.String = 'Done';
but_done.Tag = 'but_done';
but_done.Callback = @finalize;
but_done.Enable = 'off';

% initialize the gui
axes(findobj('Tag', 'ax_displ'));
imagesc(share.piv_x, share.piv_y, share.piv_m, 'AlphaData', share.piv_roi);
set(gca, 'YDir', 'Normal');
set(gca, 'Color', 0.9*[1 1 1]);
axis equal tight
colorbar

get_ctrl_pts();

hf.Visible = 'on';

end


function finalize(hui, ~)
% Finalize analysis and write results to file
% %

% make sure the last point has been recorded
next_pt([], [], 1, 0);

% get shared data, then destroy the GUI
share = get(gcf, 'UserData');

% close the GUI
close(hui.Parent);

% interpolate PIV displacements at control points
ctrl_u_piv = interp2(share.piv_x, share.piv_y, share.piv_u, share.ctrl_x, share.ctrl_y, 'linear');
ctrl_v_piv = interp2(share.piv_x, share.piv_y, share.piv_v, share.ctrl_x, share.ctrl_y, 'linear');

% compute distance to sand boundary for control points
roi = ~isnan(share.piv_u);
piv_d = bwdist(bwperim(roi));
piv_d(~roi) = NaN;
ctrl_d = interp2(share.piv_x, share.piv_y, piv_d, share.ctrl_x, share.ctrl_y, 'linear'); % pixels
pixel_to_meter = abs(share.piv_x(1)-share.piv_x(2));
ctrl_d = ctrl_d*pixel_to_meter;

% compute results vars
ctrl_m = sqrt(share.ctrl_u.^2+share.ctrl_v.^2);
ctrl_m_piv = sqrt(ctrl_u_piv.^2+ctrl_v_piv.^2);
ctrl_theta = atand(share.ctrl_v./share.ctrl_u);
ctrl_theta_piv = atand(ctrl_v_piv./ctrl_u_piv);
diff_u = 1000*(ctrl_u_piv-share.ctrl_u);
diff_v = 1000*(ctrl_v_piv-share.ctrl_v);
diff_m = 1000*(ctrl_m_piv-ctrl_m);
diff_theta = ctrl_theta_piv-ctrl_theta;

% write analysis data to netCDF file
ncid = netcdf.create(share.result_nc, 'CLOBBER');

globatt = netcdf.getConstant('NC_GLOBAL'); 
netcdf.putAtt(ncid, globatt, 'image_file_name', share.image_file);
netcdf.putAtt(ncid, globatt, 'image_file_md5', util_md5_hash(share.image_file));
netcdf.putAtt(ncid, globatt, 'displ_file_name', share.displ_file);
netcdf.putAtt(ncid, globatt, 'displ_file_md5', util_md5_hash(share.displ_file));
netcdf.putAtt(ncid, globatt, 'step_index', share.step);

dim_pt = netcdf.defDim(ncid, 'point', share.num_pts);
dim_image_x = netcdf.defDim(ncid, 'image_x', length(share.image_x));
dim_image_y = netcdf.defDim(ncid, 'image_y', length(share.image_y));
dim_piv_x = netcdf.defDim(ncid, 'piv_x', length(share.piv_x));
dim_piv_y = netcdf.defDim(ncid, 'piv_y', length(share.piv_y));

var_pt = netcdf.defVar(ncid, 'pt', 'NC_SHORT', dim_pt);
netcdf.putAtt(ncid, var_pt, 'long_name', 'control point index');
netcdf.putAtt(ncid, var_pt, 'units', '1');

var_ctrl_x = netcdf.defVar(ncid, 'ctrl_x', 'NC_DOUBLE', dim_pt);
netcdf.putAtt(ncid, var_ctrl_x, 'long_name', 'control point x-position');
netcdf.putAtt(ncid, var_ctrl_x, 'units', 'meters');

var_ctrl_y = netcdf.defVar(ncid, 'ctrl_y', 'NC_DOUBLE', dim_pt);
netcdf.putAtt(ncid, var_ctrl_y, 'long_name', 'control point y-position');
netcdf.putAtt(ncid, var_ctrl_y, 'units', 'meters');

var_ctrl_d = netcdf.defVar(ncid, 'ctrl_d', 'NC_DOUBLE', dim_pt);
netcdf.putAtt(ncid, var_ctrl_d, 'long_name', 'control point distance to boundary');
netcdf.putAtt(ncid, var_ctrl_d, 'units', 'meters');

var_ctrl_u_man = netcdf.defVar(ncid, 'ctrl_u_man', 'NC_DOUBLE', dim_pt);
netcdf.putAtt(ncid, var_ctrl_u_man, 'long_name', 'control point displacement, x-component, manual estimate');
netcdf.putAtt(ncid, var_ctrl_u_man, 'units', 'meters');

var_ctrl_v_man = netcdf.defVar(ncid, 'ctrl_v_man', 'NC_DOUBLE', dim_pt);
netcdf.putAtt(ncid, var_ctrl_v_man, 'long_name', 'control point displacement, y-component, manual estimate');
netcdf.putAtt(ncid, var_ctrl_v_man, 'units', 'meters');

var_ctrl_u_piv = netcdf.defVar(ncid, 'ctrl_u_piv', 'NC_DOUBLE', dim_pt);
netcdf.putAtt(ncid, var_ctrl_u_piv, 'long_name', 'control point displacement, x-component, PIV estimate');
netcdf.putAtt(ncid, var_ctrl_u_piv, 'units', 'meters');

var_ctrl_v_piv = netcdf.defVar(ncid, 'ctrl_v_piv', 'NC_DOUBLE', dim_pt);
netcdf.putAtt(ncid, var_ctrl_v_piv, 'long_name', 'control point displacement, y-component, PIV estimate');
netcdf.putAtt(ncid, var_ctrl_v_piv, 'units', 'meters');

var_image_x = netcdf.defVar(ncid, 'image_x', 'NC_DOUBLE', dim_image_x);
netcdf.putAtt(ncid, var_image_x, 'long_name', 'image grid x-coordinate');
netcdf.putAtt(ncid, var_image_x, 'units', 'meters');

var_image_y = netcdf.defVar(ncid, 'image_y', 'NC_DOUBLE', dim_image_y);
netcdf.putAtt(ncid, var_image_y, 'long_name', 'image grid y-coordinate');
netcdf.putAtt(ncid, var_image_y, 'units', 'meters');

var_ini = netcdf.defVar(ncid, 'ini', 'NC_DOUBLE', [dim_image_y, dim_image_x]);
netcdf.putAtt(ncid, var_ini, 'long_name', 'image at initial time');
netcdf.putAtt(ncid, var_ini, 'units', 'normalized intensity');

var_fin = netcdf.defVar(ncid, 'fin', 'NC_DOUBLE', [dim_image_y, dim_image_x]);
netcdf.putAtt(ncid, var_fin, 'long_name', 'image at final time'); 
netcdf.putAtt(ncid, var_fin, 'units', 'normalized intensity');

var_piv_x = netcdf.defVar(ncid, 'piv_x', 'NC_DOUBLE', dim_piv_x);
netcdf.putAtt(ncid, var_piv_x, 'long_name', 'PIV grid x-coordinate');
netcdf.putAtt(ncid, var_piv_x, 'units', 'meters');

var_piv_y = netcdf.defVar(ncid, 'piv_y', 'NC_DOUBLE', dim_piv_y);
netcdf.putAtt(ncid, var_piv_y, 'long_name', 'PIV grid y-coordinate');
netcdf.putAtt(ncid, var_piv_y, 'units', 'meters');

var_piv_u = netcdf.defVar(ncid, 'piv_u', 'NC_DOUBLE', [dim_piv_y, dim_piv_x]);
netcdf.putAtt(ncid, var_piv_u, 'long_name', 'displacement, x-component, piv estimate');
netcdf.putAtt(ncid, var_piv_u, 'units', 'meters/step');

var_piv_v = netcdf.defVar(ncid, 'piv_v', 'NC_DOUBLE', [dim_piv_y, dim_piv_x]);
netcdf.putAtt(ncid, var_piv_v, 'long_name', 'displacement, y-component, piv estimate');
netcdf.putAtt(ncid, var_piv_v, 'units', 'meters/step');

netcdf.endDef(ncid);

netcdf.putVar(ncid, var_pt, 1:share.num_pts);
netcdf.putVar(ncid, var_ctrl_x, share.ctrl_x);
netcdf.putVar(ncid, var_ctrl_y, share.ctrl_y);
netcdf.putVar(ncid, var_ctrl_d, ctrl_d);
netcdf.putVar(ncid, var_ctrl_u_man, share.ctrl_u);
netcdf.putVar(ncid, var_ctrl_v_man, share.ctrl_v);
netcdf.putVar(ncid, var_ctrl_u_piv, ctrl_u_piv);
netcdf.putVar(ncid, var_ctrl_v_piv, ctrl_v_piv);
netcdf.putVar(ncid, var_image_x, share.image_x);
netcdf.putVar(ncid, var_image_y, share.image_y);
netcdf.putVar(ncid, var_ini, share.ini);
netcdf.putVar(ncid, var_fin, share.fin);
netcdf.putVar(ncid, var_piv_x, share.piv_x);
netcdf.putVar(ncid, var_piv_y, share.piv_y);
netcdf.putVar(ncid, var_piv_u, share.piv_u);
netcdf.putVar(ncid, var_piv_v, share.piv_v);

netcdf.close(ncid);

% plot control point location on displacement magnitude, save
hf = figure;
hf.Units = 'Normalized';
hf.Position = [0, 0.3, 1, 0.3];
hf.Name = 'piv_check_manual';

ax = axes();
imagesc(share.piv_x, share.piv_y, share.piv_m*1000, 'AlphaData', ~isnan(share.piv_u));
hcb = colorbar;
hcb.Label.String = 'mm/step';
ax.YDir = 'Normal';
hold on
plot(share.ctrl_x, share.ctrl_y, '*k');
xlabel('x-position, m');
ylabel('y-position, m');
title('Location of control points for manual PIV check');

saveas(hf, share.result_fig_loc);

% plot PDF of piv-manual
hf = figure;
hf.Name = 'piv_check_manual';

subplot(2,2,1);
[y, x] = ksdensity(diff_u);
plot(x, y);
xlabel('u_{PIV-manual}, mm');
ylabel('\rho');

subplot(2,2,2);
[y, x] = ksdensity(diff_v);
plot(x, y);
xlabel('v_{PIV-manual}, mm');
ylabel('\rho');

subplot(2,2,3);
[y, x] = ksdensity(diff_m);
plot(x, y);
xlabel('mag_{PIV-manual}, mm');
ylabel('\rho');

subplot(2,2,4);
[y, x] = ksdensity(diff_theta);
plot(x, y);
xlabel('\theta_{PIV-manual}, deg');
ylabel('\rho');

saveas(hf, share.result_fig_pdf);

% plot piv vs. manual (u, v, mag, theta)
hf = figure;
hf.Name = 'piv_check_manual';

subplot(2,2,1);
plot(1000*ctrl_d, diff_u, 'xk');
xlabel('Dist to Bnd, mm');
ylabel('u_{PIV-manual}, mm');

subplot(2,2,2);
plot(1000*ctrl_d, diff_v, 'xk');
xlabel('Dist to Bnd, mm');
ylabel('v_{PIV-manual}, mm');

subplot(2,2,3);
plot(1000*ctrl_d, diff_m, 'xk');
xlabel('Dist to Bnd, mm');
ylabel('mag_{PIV-manual}, mm');

subplot(2,2,4);
plot(1000*ctrl_d, diff_theta, 'xk');
xlabel('Dist to Bnd, mm');
ylabel('\theta_{PIV-manual}, deg');

saveas(hf, share.result_fig_vs);

end

function delete_pt(~, ~)
% Delete the current point

% get shared data
hfig = gcf();
share = get(gcf, 'UserData');

% delete point
share.num_pts = share.num_pts-1;
share.ctrl_x(share.ii) = [];
share.ctrl_y(share.ii) = [];
share.ctrl_u(share.ii) = [];
share.ctrl_v(share.ii) = [];

% set shared data
set(hfig, 'UserData', share);

% refresh the GUI for the next point
next_pt([], [], 0, 0);

end

function [share] = record_pt(share)
% Record analysis results from current axes (u_man, v_man, new ctrl_x, ctrl_y)
% %

% get location of initial and final impoints
axes(findobj('Tag', 'ax_ini'));
hpt = get(gca, 'UserData');
ini_xy = hpt.getPosition();

axes(findobj('Tag', 'ax_fin'));
hpt = get(gca, 'UserData');
fin_xy = hpt.getPosition();

% compute manually estimated displacement
man_uv = fin_xy-ini_xy;
share.ctrl_u(share.ii) = man_uv(1);
share.ctrl_v(share.ii) = man_uv(2);

% update location of the control point to the midpoint
mid_xy = mean([ini_xy; fin_xy]);
share.ctrl_x(share.ii) = mid_xy(1);
share.ctrl_y(share.ii) = mid_xy(2);

end


function next_pt(~, ~, record, step)
% Record results from the last point (record == 1) or don't (record == 0), and
% analyze the next (step == 1) or previous (step == -1) control point
% %

% constants
dim = [-0.003, 0.003];

% get shared data
hfig = gcf();
share = get(gcf, 'UserData');


% record results from last point
if record == 1
    disp([share.num_pts, share.ii, share.ctrl_x(share.ii), share.ctrl_y(share.ii)]);
    share = record_pt(share);
    disp([share.num_pts, share.ii, share.ctrl_x(share.ii), share.ctrl_y(share.ii)]);
end

% change to next/prev control point, wrapping around as needed
share.ii = share.ii+step;
if share.ii > share.num_pts; share.ii = 1;             end
if share.ii < 1;             share.ii = share.num_pts; end

% update counter
counter = findobj('Tag', 'text_curr_pt');
counter.String = sprintf('Point %i of %i', share.ii, share.num_pts); 

% update ini and fin plots
axes(findobj('Tag', 'ax_ini'));
hpt = get(gca, 'UserData');
guess_x = share.ctrl_x(share.ii) - 0.5*share.ctrl_u(share.ii); % use manual estimate so I can refine
guess_y = share.ctrl_y(share.ii) - 0.5*share.ctrl_v(share.ii);

hpt.setPosition([guess_x, guess_y]);
xlim(guess_x+dim);
ylim(guess_y+dim);

axes(findobj('Tag', 'ax_fin'));
hpt = get(gca, 'UserData');
guess_x = share.ctrl_x(share.ii) + 0.5*share.ctrl_u(share.ii);
guess_y = share.ctrl_y(share.ii) + 0.5*share.ctrl_v(share.ii);
hpt.setPosition([guess_x, guess_y]);
xlim(guess_x+dim);
ylim(guess_y+dim);

% set shared data
set(hfig, 'UserData', share);

end

function start_analysis(hui, ~)
% Switch from initialization to analysis mode
% %

% disable / enable controls
hui.Enable = 'off';
h = findobj('Tag', 'edit_num_pts'); h.Enable = 'off';
h = findobj('Tag', 'but_get_pts');  h.Enable = 'off';
h = findobj('Tag', 'but_next_pt'); h.Enable = 'on';
h = findobj('Tag', 'but_prev_pt'); h.Enable = 'on';
h = findobj('Tag', 'but_delete'); h.Enable = 'on';
h = findobj('Tag', 'but_done'); h.Enable = 'on';

% start analysis for first point
next_pt([], [], 0, 1);

end

function get_ctrl_pts(~, ~)
% Get randomly distributed control points
% %

% get shared data
share = get(gcf, 'UserData');

% (re)allocate control point vectors 
share.ctrl_x = nan(share.num_pts, 1);
share.ctrl_y = nan(share.num_pts, 1);
share.ctrl_u = nan(share.num_pts, 1);
share.ctrl_v = nan(share.num_pts, 1);

% generate num_pts random points within the data roi
get_pt_x = @() rand(1)*range(share.piv_x)+min(share.piv_x);
get_pt_y = @() rand(1)*range(share.piv_y)+min(share.piv_y);

for ii = 1:share.num_pts
    while isnan(share.ctrl_x(ii))
       pt_x = get_pt_x();
       pt_y = get_pt_y();
       if interp2(share.piv_x, share.piv_y, share.piv_roi, pt_x, pt_y, 'linear') == 1
           share.ctrl_x(ii) = pt_x;
           share.ctrl_y(ii) = pt_y;
           break
       end
    end
end

% interpolate piv displacements at the control points, initial guess from PIV results
share.ctrl_u = interp2(share.piv_x, share.piv_y, share.piv_u, share.ctrl_x, share.ctrl_y, 'linear');
share.ctrl_v = interp2(share.piv_x, share.piv_y, share.piv_v, share.ctrl_x, share.ctrl_y, 'linear');

% TO DO: compute distance to the boundary at the control points

% set shared data
set(gcf, 'UserData', share);

% update displacement plot
update_ax_displ;

end

function update_ax_displ()
% Replot the location of control points in the displacement magnitude plot
% %

hax = findobj('Tag', 'ax_displ');
axes(hax);
delete(findobj(hax, 'Type', 'Line')); % remove any existing points
hold on
plot(hax.Parent.UserData.ctrl_x, hax.Parent.UserData.ctrl_y, 'xk');
hold off

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



