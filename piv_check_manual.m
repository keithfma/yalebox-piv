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

% sanity check
validateattributes(image_file,  {'char'}, {'vector'}, mfilename, 'image_file');
validateattributes(displ_file,  {'char'}, {'vector'}, mfilename, 'displ_file');
validateattributes(result_file, {'char'}, {'vector'}, mfilename, 'result_file');
validateattributes(step, {'numeric'}, {'scalar'}, mfilename, 'step');

% create struct for shared data
share = struct(...
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
    'ctrl_u_piv', [], ...
    'ctrl_v_piv', [], ...
    'ctrl_u_man', [], ...
    'ctrl_v_man', [], ...
    'ctrl_d', [], ...
    'ii', 0);

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
set(gca, 'Tag', 'ax_ini', 'NextPlot','add', 'YDir', 'Normal', 'UserData', hpt);
title('Initial Image');

axes('Position', [0.45, 0.05, 0.35, 0.5]);
imagesc(share.piv_x, share.piv_y, share.fin);
hpt = impoint(gca, [0, 0]);
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
but_next_pt.Callback = {@next_pt, 1};
but_next_pt.Enable = 'off';

but_prev_pt = uicontrol('Style', 'pushbutton');
but_prev_pt.Units = 'Normalized';
but_prev_pt.Position = [0.85, 0.35, 0.1, 0.05];
but_prev_pt.String = 'Previous Point';
but_prev_pt.Tag = 'but_prev_pt';
but_prev_pt.Callback = {@next_pt, -1};
but_prev_pt.Enable = 'off';

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

% function record_pt()
% Record analysis results from current axes (u_man, v_man, new ctrl_x, ctrl_y)

function next_pt(~, ~, step)
% Record results from the last point, and analyze the next (step == 1) or
% previous (step == -1) control point
% %

% get shared data
hfig = gcf();
share = get(gcf, 'UserData');

% record results from last point
if share.ii ~= 0 % skip first time
    % call a "record_results" function...
end

% change to next/prev control point, wrapping around as needed
share.ii = share.ii+step;
if share.ii > share.num_pts; share.ii = 1; end
if share.ii < 1; share.ii = share.num_pts; end

% update ini and fin plots
axes(findobj('Tag', 'ax_ini'));
hpt = get(gca, 'UserData');
guess_x = share.ctrl_x(share.ii) - 0.5*share.ctrl_u_piv(share.ii);
guess_y = share.ctrl_y(share.ii) - 0.5*share.ctrl_v_piv(share.ii);
hpt.setPosition([guess_x, guess_y]);
xlim(guess_x+[-0.005, 0.005]);
ylim(guess_y+[-0.005, 0.005]);

axes(findobj('Tag', 'ax_fin'));
hpt = get(gca, 'UserData');
guess_x = share.ctrl_x(share.ii) + 0.5*share.ctrl_u_piv(share.ii);
guess_y = share.ctrl_y(share.ii) + 0.5*share.ctrl_v_piv(share.ii);
hpt.setPosition([guess_x, guess_y]);
xlim(guess_x+[-0.005, 0.005]);
ylim(guess_y+[-0.005, 0.005]);

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

% start analysis for first point
next_pt([], [], 1);

end

function get_ctrl_pts(~, ~)
% Get randomly distributed control points
% %

% get shared data
share = get(gcf, 'UserData');

% (re)allocate control point vectors 
share.ctrl_x     = nan(share.num_pts, 1);
share.ctrl_y     = nan(share.num_pts, 1);
share.ctrl_u_piv = nan(share.num_pts, 1);
share.ctrl_v_piv = nan(share.num_pts, 1);
share.ctrl_u_man = nan(share.num_pts, 1);
share.ctrl_v_man = nan(share.num_pts, 1);
share.ctrl_d     = nan(share.num_pts, 1);

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

% interpolate piv displacements at the control points
share.ctrl_u_piv = interp2(share.piv_x, share.piv_y, share.piv_u, share.ctrl_x, share.ctrl_y, 'linear');
share.ctrl_v_piv = interp2(share.piv_x, share.piv_y, share.piv_v, share.ctrl_x, share.ctrl_y, 'linear');

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



