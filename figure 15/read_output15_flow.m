% Read Nek5000 output data
clc; clear all;close all;

folderpath = './'

u_mean = readmatrix("u_mean_nek.csv");
v_mean = readmatrix("v_mean_nek.csv");
x_mean = readmatrix("x_mean_nek.csv")
y_mean = readmatrix("y_mean_nek.csv")
nlevels = 40;

rwb1 = bluewhitered;

%% U mean nek

fig = figure;
set(fig, 'Position', [100, 100, 1000, 800]);
contourf(x_mean, y_mean, u_mean, nlevels);
shading interp
colormap(rwb1);

% -------- Colorbar setup (unchanged, but cleaned) ----------
c = colorbar;
%c.Limits = [min(u_2D_mean(:)) max(u_2D_mean(:))];
c.Limits = [-0.24 1.98];
numTicks = 6;
minVal = c.Limits(1);
maxVal = c.Limits(2);
ticks   = linspace(minVal, maxVal, numTicks);

c.Ticks      = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
c.FontSize   = 32;

fprintf('U Velocity Nek min: %.3f\n', minVal);
fprintf('U Velocity Nek max: %.3f\n', maxVal);

set(gca, 'FontSize', 32);
xlabel('x', 'FontSize', 40);
ylabel('y', 'FontSize', 40);

% ================= AXIS TICKS FROM 0 ON BOTH AXES =================
ax = gca;

% ---- X-axis: from 0 to max, like snapshots code ----
xmin = min(x_mean(:));
xmax = max(x_mean(:));

% if your domain is [0,9], enforce that explicitly:
ax.XLim      = [0 9];       % or [xmin xmax] if you prefer data-based
ax.XTickMode = 'manual';
ax.XTick     = 0:3:9;       % ticks at 0,3,6,9
ax.XTickLabel = arrayfun(@(x) sprintf('%d', x), ax.XTick, 'UniformOutput', false);

% ---- Y-axis: from 0 to 3, independent from x-axis ----
ymin = min(y_mean(:));
ymax = max(y_mean(:));

% if your physical wall-normal extent is [0,3], enforce:
ax.YLim      = [0 3];       % or [ymin ymax] if you'd rather auto
ax.YTickMode = 'manual';
ax.YTick     = 0:1:3;       % ticks at 0,1,2,3
ax.YTickLabel = arrayfun(@(y) sprintf('%d', y), ax.YTick, 'UniformOutput', false);
% =====================================================

saveas(gcf, 'U_nek_Re190.png');

%% V nek ------------------------------------------------------------------
fig = figure;
set(fig, 'Position', [100, 100, 1000, 800]);

contourf(x_mean, y_mean, v_mean, nlevels);
shading interp
colormap(rwb1);

% ---------------- COLORBAR ----------------
c = colorbar;
%c.Limits = [min(v_2D_mean(:)) max(v_2D_mean(:))];
c.Limits = [-0.19 0.22];
numTicks = 6;

minVal = c.Limits(1);
maxVal = c.Limits(2);
ticks  = linspace(minVal, maxVal, numTicks);

c.Ticks      = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
c.FontSize   = 32;

fprintf('V Velocity Colorbar Min: %.3f\n', minVal);
fprintf('V Velocity Colorbar Max: %.3f\n', maxVal);

% ---------------- AXES ----------------
ax = gca;
set(ax, 'FontSize', 32);

xlabel('x', 'FontSize', 40);
ylabel('y', 'FontSize', 40);

% ================= AXIS TICKS FROM 0 ON BOTH AXES =================
% ---- X-axis ----
ax.XLim      = [0 9];
ax.XTickMode = 'manual';
ax.XTick     = 0:3:9;
ax.XTickLabel = arrayfun(@(x) sprintf('%d', x), ax.XTick, 'UniformOutput', false);

% ---- Y-axis ----
ax.YLim      = [0 3];
ax.YTickMode = 'manual';
ax.YTick     = 0:1:3;
ax.YTickLabel = arrayfun(@(y) sprintf('%d', y), ax.YTick, 'UniformOutput', false);
% =================================================================

saveas(gcf, 'V_nek_Re190.png');
%% ------------------------------------------------------------------------