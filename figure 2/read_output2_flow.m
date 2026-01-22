clc; clear all;close all;

folderpath = './'

X = readmatrix('X.csv');
Y = readmatrix('Y.csv');

x_solid = readmatrix('x_solid.csv');
y_solid = readmatrix('y_solid.csv');
Lx = 9.0;
epsilon = 0.5;
Ly = 1.0 + epsilon;

u_mean_zt = readmatrix('dedalus_u_mean_zt.csv');
v_mean_zt = readmatrix('dedalus_v_mean_zt.csv');

tick_enforcing = "manual";
nLevels = 30;
rwb1 = bluewhitered
%% =========================
%   u_mean_zt contour (single)
% =========================
% rwb1 = bluewhitered;
fig_uMean = figure(101);
clf(fig_uMean)
set(fig_uMean, 'Position', [100, 100, 1000, 800]);

contourf(X, Y, u_mean_zt, nLevels);
colormap(rwb1);

% === Y-AXIS RELABELING (+1.5 shift) ===
yt = get(gca, 'YTick');                    % original ticks (e.g., -1.5..1.5)
set(gca, 'YTickLabel', yt + 1.5);          % relabeled as 0..3

c = colorbar;
c.FontSize = 32;

% Fixed number of ticks and formatted labels
numTicks = 6;
minVal = min(c.Limits);
maxVal = max(c.Limits);
ticks = linspace(minVal, maxVal, numTicks);
filteredTicks = ticks(ticks ~= minVal & ticks ~= maxVal);
%c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), filteredTicks, 'UniformOutput', false);

c.Ticks = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
c.FontSize = 32;

fprintf('U Velocity Snapshot %d Min: %.3f\n', i, minVal);
fprintf('U Velocity Snapshot %d Max: %.3f\n', i, maxVal);

%title(sprintf("U Velocity (contour) - Snapshot %d (t = %d)", i, snap_idx(i)), 'FontSize', 22)
set(gca, 'FontSize', 32);
xlabel('x', 'FontSize', 40);
ylabel('y', 'FontSize', 36);

hold on
plot(x_solid, y_solid, '.', 'Color', 'white', 'MarkerSize', 12);
hold off

% ===========================================================
%      ADD MIN & MAX TICKS TO BOTH X AND Y AXES HERE
% ===========================================================
ax = gca;

% X-axis
xmin = ax.XLim(1);
xmax = ax.XLim(2);
origX = ax.XTick;

if tick_enforcing == "automatic"
    newX = unique([xmin; origX(:); xmax]);
else
    ax.XLim = [0 9];              % force exact endpoints
    ax.XTickMode = 'manual';      % prevent auto changes
    newX = 0:3:9;
end

ax.XTick = newX;
ax.XTickLabel = arrayfun(@(x) sprintf('%d', x), newX, 'UniformOutput', false);

% Y-axis (before +1.5 shift, so use actual coords)
ymin = ax.YLim(1);
ymax = ax.YLim(2);
origY = ax.YTick;

if tick_enforcing == "automatic"
    newY = unique([ymin; origY(:); ymax]);
else
    ax.YLim = [-1.5 1.5];              % force exact endpoints
    ax.YTickMode = 'manual';      % prevent auto changes
    newY = -1.5:1:1.5;
end

ax.YTick = newY;
% apply the +1.5 shift to labels
ax.YTickLabel = arrayfun(@(y) sprintf('%d', y + 1.5), newY, 'UniformOutput', false);
% ===========================================================

saveas(fig_uMean,'u_mean_zt_contour.png');



%% =========================
%   v_mean_zt contour (single)
% =========================
fig_vMean = figure(102);
clf(fig_vMean)
set(fig_vMean, 'Position', [100, 100, 1000, 800]);

contourf(X, Y, v_mean_zt, nLevels);
colormap(rwb1);

% === Y-AXIS RELABELING (+1.5 shift) ===
yt = get(gca, 'YTick');                    % original ticks (e.g., -1.5..1.5)
set(gca, 'YTickLabel', yt + 1.5);          % relabeled as 0..3

c = colorbar;
c.FontSize = 32;

% Fixed number of ticks and formatted labels
numTicks = 6;
minVal = min(c.Limits);
maxVal = max(c.Limits);
ticks = linspace(minVal, maxVal, numTicks);
filteredTicks = ticks(ticks ~= minVal & ticks ~= maxVal);
%     c.Ticks = filteredTicks;
%     c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), filteredTicks, 'UniformOutput', false);
c.Ticks = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
c.FontSize = 32;

fprintf('V Velocity Snapshot %d Min: %.3f\n', i, minVal);
fprintf('V Velocity Snapshot %d Max: %.3f\n', i, maxVal);

%title(sprintf("V Velocity (contour) - Snapshot %d (t = %d)", i, snap_idx(i)), 'FontSize', 22)
set(gca, 'FontSize', 32);
xlabel('x', 'FontSize', 40);
ylabel('y', 'FontSize', 40);

hold on
plot(x_solid, y_solid, '.', 'Color', 'white', 'MarkerSize', 12);
hold off

% ===========================================================
%      ADD MIN & MAX TICKS TO BOTH X AND Y AXES HERE
% ===========================================================
ax = gca;

% X-axis
xmin = ax.XLim(1);
xmax = ax.XLim(2);
origX = ax.XTick;

if tick_enforcing == "automatic"
    newX = unique([xmin; origX(:); xmax]);
else
    ax.XLim = [0 Lx];              % force exact endpoints
    ax.XTickMode = 'manual';      % prevent auto changes
    newX = 0:3:9;
end

ax.XTick = newX;
ax.XTickLabel = arrayfun(@(x) sprintf('%d', x), newX, 'UniformOutput', false);

% Y-axis (before +1.5 shift, so use actual coords)
ymin = ax.YLim(1);
ymax = ax.YLim(2);
origY = ax.YTick;

if tick_enforcing == "automatic"
    newY = unique([ymin; origY(:); ymax]);
else
    ax.YLim = [-1.5 1.5];              % force exact endpoints
    ax.YTickMode = 'manual';      % prevent auto changes
    newY = -1.5:1:1.5;
end

ax.YTick = newY;
% apply the +1.5 shift to labels
ax.YTickLabel = arrayfun(@(y) sprintf('%d', y + 1.5), newY, 'UniformOutput', false);
% ===========================================================

saveas(fig_vMean,'v_mean_zt_contour.png');

writematrix(u_mean_zt,'u_mean_zt.csv');
writematrix(v_mean_zt,'v_mean_zt.csv');