function display_image_pair(xx, yy, ini, ini_roi, fin, fin_roi, which, fill)
% function display_image_pair(xx, yy, ini, ini_roi, fin, fin_roi, which, fill)
%
% Simple plot of image pair with area outside ROI set to 0
%
% Arguments:
%   xx, yy: Vector, world coordinate vectors for image matrices
%   ini, fin: 2D matrix, intensity values for initial and final images 
%   ini_roi, fin_roi: 2/3D matrix, logical flags for inital and final images,
%       indicating if a pixel is sand (1) or not (0)
%   which: Optional string, select which kind of display, options are
%       'subplot' and 'toggle' 
%   fill: Optional, set true to fill non-sand pixels as in PIV routine
%   color: Optional, set false to display equalized grayscale, else will
%       show rgb color (default)
% %

% set defaults
if nargin < 7; which = 'subplot'; end
if nargin < 8; fill = false; end

if fill
    % use grayscale, and apply padding
    update_path('piv')
    pad_width = 25; % constant pad width because it does not seem worth specifying
    xx0 = xx; % coordinates get overwritten, use the originals
    yy0 = yy;
    [xx, yy, ini] = piv_fill(xx0, yy0, ini, ini_roi, pad_width, pad_width);
    [~, ~, fin] = piv_fill(xx0, yy0, fin, fin_roi, pad_width, pad_width);

else
    % apply mask
    if size(ini, 3) == 3; ini_roi = repmat(ini_roi, 1, 1, 3); end
    ini(~ini_roi) = 0;
    
    if size(fin, 3) == 3; fin_roi = repmat(fin_roi, 1, 1, 3); end
    fin(~fin_roi) = 0;

end

% setup figure for specified display type
hf = figure;

if strcmp(which, 'subplot')
    ax_ini = subplot(2,1,1);
    ax_fin = subplot(2,1,2);
    linkaxes([ax_ini, ax_fin]);
    
elseif strcmp(which, 'toggle')
    ax_ini = axes();
    ax_fin = axes();
    
    switch_btn = uicontrol(...
        'Style', 'pushbutton', ...
        'String', 'Switch Image',...
        'Callback', @switch_axes, ...
        'Units', 'Normalized');
    switch_btn.Position = [0.4, switch_btn.Position(2), ...
                           0.2, switch_btn.Position(4)];
    
else
    error('Invalid argument value: which = %s', which);
end

% plot images
axes(ax_ini);
imagesc(xx, yy, ini);
colormap('gray');
ax_ini.YDir = 'normal';
axis equal tight
title('Initial Image');

axes(ax_fin);
imagesc(xx, yy, fin);
colormap('gray');
ax_fin.YDir = 'normal';
axis equal tight
title('Final Image');

linkaxes([ax_ini, ax_fin]);

% some final formatting
addToolbarExplorationButtons(gcf);
ax_ini.Toolbar.Visible = 'off';
ax_fin.Toolbar.Visible = 'off';

if strcmp(which, 'subplot')
    hf.Units = 'Normalized';
    hf.OuterPosition = [0, 0, 1, 1];

elseif strcmp(which, 'toggle')
    hf.Units = 'Normalized';
    hf.OuterPosition = [0, hf.OuterPosition(2), 1, hf.OuterPosition(4)];
    toggle_axis(ax_fin);

end

end

function toggle_axis(this_ax)
% make specific axis and its children visible/invisible

if strcmp(this_ax.Visible, 'on')
    change_to = 'off';
    uistack(this_ax, 'top');
elseif strcmp(this_ax.Visible, 'off')
    change_to = 'on';
else
    error('how did this happen, anyway?');
end

for jj = 1:length(this_ax.Children)
    this_ax.Children(jj).Visible = change_to;
end
this_ax.Visible = change_to;

end

function switch_axes(~, ~)
% callback that switches axes

ax_list = findobj(gcf, 'type', 'axes');

for ii = 1:length(ax_list)
    toggle_axis(ax_list(ii));
end

end
