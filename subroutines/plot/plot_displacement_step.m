function plot_displacement_step(piv_file, step)
% function plot_displacement_step(piv_file, step)
%
% Display displacements for specified step from PIV data file
%
% Arguments:
%   piv_file: string, path to PIV data file, as created by piv.m
%   step: int, step (not step index) to be displayed
% %

% get relevant data
data = matfile(piv_file, 'Writable', false);
idx = find(data.step == step);
x = data.x*100;  % m -> cm
y = data.y*100;  % m -> cm
u = data.u(:, :, idx)*1000;  % m -> mm
v = data.v(:, :, idx)*1000;  % m -> mm
m = sqrt(u.^2 + v.^2);
roi = ~isnan(u) & ~isnan(v);

% get data range
[data_i, ~] = find(roi);
ylim = [min(y), max(y(data_i)) + 0.05*range(y(data_i))];
xlim = [min(x), max(x)];

% define formatting function
function format_subplot(hax, him, val, ttl, units)
    him.AlphaData = roi;
    caxis(hax, prctile(val(roi), [5, 95]));
    hax.YDir = 'Normal';
    hcb = colorbar('EastOutside');
    hcb.XLabel.String = units;
    axis equal
    hax.XLim = xlim;
    hax.YLim = ylim;
    hax.FontSize = 12;
    hax.XLabel.String = 'x [cm]';
    hax.YLabel.String = 'y [cm]';
    hax.Title.String = ttl;
end

% create figure, fullscreen
hf = figure('units','normalized','outerposition',[0 0 1 1]);
hf.Color = 'w';
hf.Name = sprintf('PIV Displacements - Step %d', step);

% subplot: displacement magnitude with (downsampled) vector arrows
hax_m = subplot(3, 1, 1);
him_m = imagesc([x(1), x(end)], [y(1), y(end)], m);
hold on 
q = checkermask(m, 5);
[xg, yg] = meshgrid(x, y);
quiver(xg(q), yg(q), u(q), v(q), 'Color', 'k');
format_subplot(hax_m, him_m, m, 'Displacement Magnitide and Direction', 'mm/step');

% subplot: horizontal displacement magnitude
hax_u = subplot(3, 1, 2);
him_u = imagesc([x(1), x(end)], [y(1), y(end)], u);
format_subplot(hax_u, him_u, u, 'Horizontal Displacement Magnitide', 'mm/step');

% subplot: vertical displacement magnitude
hax_v = subplot(3, 1, 3);
him_v = imagesc([x(1), x(end)], [y(1), y(end)], v);
format_subplot(hax_v, him_v, v, 'Vertical Displacement Magnitide', 'mm/step');

fprintf('DONE\n');

end

