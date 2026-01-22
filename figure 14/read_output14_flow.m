% Read Nek5000 output data
clc; clear all;close all;

folderpath = './'

u_2D_mean = readmatrix("u_mean_nek_Re100.csv");
v_2D_mean = readmatrix("v_mean_nek_Re100.csv");
x_2D = readmatrix("x_mean_nek_Re100.csv")
y_2D = readmatrix("y_mean_nek_Re100.csv")

rwb1 = bluewhitered;
nlevels = 40;
%% U nek
fig = figure;
%set(fig, 'Position', [100, 100, 900, 750]);  % Resize figure
set(fig, 'Position', [100, 100, 1000, 800]);
contourf(x_2D, y_2D, u_2D_mean, nlevels);
shading interp
colormap(rwb1);

c = colorbar;
existingTicks = c.Ticks;
numTicks = 6;
c.Limits = [-0.15 1.93]
minVal = min(c.Limits);
maxVal = max(c.Limits);
ticks = linspace(minVal, maxVal, numTicks);% Remove min and max from ticks
filteredTicks = existingTicks(existingTicks ~= minVal & existingTicks ~= maxVal);
% c.Ticks = filteredTicks;
%c.Ticks = linspace(c.Limits(1), c.Limits(2), numTicks);
%c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), c.Ticks, 'UniformOutput', false);
c.Ticks = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
c.FontSize = 32;

% Print min/max separately
fprintf('U Velocity Nek min: %.3f\n', minVal);
fprintf('V Velocity Nek max: %.3f\n', maxVal);

set(gca, 'FontSize', 32);
xlabel('X', 'FontSize', 32);
ylabel('Y', 'FontSize', 32);
%title("U Velocity nek5000", 'FontSize', 18);

hold on
ax = gca;

xmin = min(x_2D(:));
xmax = max(x_2D(:));

current = ax.XTick;
newticks = unique([xmin; current(:); xmax]);

ax.XTick = newticks;
ax.XTickLabel = arrayfun(@(x) sprintf('%.1f', x), newticks, 'UniformOutput', false);
saveas(gcf, 'U_nek_Re100.png');


%% V nek ------------------------------------------------------------------
fig = figure;
%set(fig, 'Position', [100, 100, 900, 750]);  % Resize figure
set(fig, 'Position', [100, 100, 1000, 800]);
contourf(x_2D, y_2D, v_2D_mean, nlevels);
shading interp
colormap(rwb1);

c = colorbar;
existingTicks = c.Ticks;
numTicks = 6;
c.Limits = [-0.07 0.28]
minVal = min(c.Limits);
maxVal = max(c.Limits);
ticks = linspace(minVal, maxVal, numTicks);% Remove min and max from ticks
% Remove min and max from ticks
filteredTicks = existingTicks(existingTicks ~= minVal & existingTicks ~= maxVal);

c.Ticks = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
c.FontSize = 32;

% Print min/max separately
fprintf('V Velocity Colorbar Min: %.3f\n', minVal);
fprintf('V Velocity Colorbar Max: %.3f\n', maxVal);

set(gca, 'FontSize', 32);
xlabel('X', 'FontSize', 32);
ylabel('Y', 'FontSize', 32);
%title("V Velocity nek5000", 'FontSize', 18);

hold on
ax = gca;

xmin = min(x_2D(:));
xmax = max(x_2D(:));

current = ax.XTick;
newticks = unique([xmin; current(:); xmax]);

ax.XTick = newticks;
ax.XTickLabel = arrayfun(@(x) sprintf('%.1f', x), newticks, 'UniformOutput', false);
saveas(gcf, 'V_nek_Re100.png');