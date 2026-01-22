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
yt = get(gca, 'YTick');
set(gca, 'YTickLabel', yt + 1.5);

c = colorbar;
c.FontSize = 32;

numTicks = 6;
minVal = min(c.Limits);
maxVal = max(c.Limits);
ticks = linspace(minVal, maxVal, numTicks);

c.Ticks = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
c.FontSize = 32;

fprintf('U_MEAN Min: %.3f\n', minVal);
fprintf('U_MEAN Max: %.3f\n', maxVal);

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
    ax.XLim = [0 9];
    ax.XTickMode = 'manual';
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
    ax.YLim = [-1.5 1.5];
    ax.YTickMode = 'manual';
    newY = -1.5:1:1.5;
end
ax.YTick = newY;
ax.YTickLabel = arrayfun(@(y) sprintf('%d', y + 1.5), newY, 'UniformOutput', false);
% ===========================================================

saveas(fig_uMean, fullfile(contourf_folder, 'u_mean_zt_contour.png'));


%% =========================
%   v_mean_zt contour (single)
% =========================
fig_vMean = figure(102);
clf(fig_vMean)
set(fig_vMean, 'Position', [100, 100, 1000, 800]);

contourf(X, Y, v_mean_zt, nLevelsv);
colormap(rwb1);

% === Y-AXIS RELABELING (+1.5 shift) ===
yt = get(gca, 'YTick');
set(gca, 'YTickLabel', yt + 1.5);

c = colorbar;
c.FontSize = 32;

numTicks = 6;
minVal = min(c.Limits);
maxVal = max(c.Limits);
ticks = linspace(minVal, maxVal, numTicks);

c.Ticks = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
c.FontSize = 32;

fprintf('V_MEAN Min: %.3f\n', minVal);
fprintf('V_MEAN Max: %.3f\n', maxVal);

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
    ax.XLim = [0 Lx];
    ax.XTickMode = 'manual';
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
    ax.YLim = [-1.5 1.5];
    ax.YTickMode = 'manual';
    newY = -1.5:1:1.5;
end
ax.YTick = newY;
ax.YTickLabel = arrayfun(@(y) sprintf('%d', y + 1.5), newY, 'UniformOutput', false);
% ===========================================================

saveas(fig_vMean, fullfile(contourf_folder, 'v_mean_zt_contour.png'));


%% =========================
%   magnitude_mean contour (single)
% =========================
% rwb1 = bluewhitered;
fig_magMean = figure(103);
clf(fig_magMean)
set(fig_magMean, 'Position', [100, 100, 1000, 800]);

contourf(X, Y, magnitude_mean, nLevels);
colormap(rwb1);

% === Y-AXIS RELABELING (+1.5 shift) ===
yt = get(gca, 'YTick');
set(gca, 'YTickLabel', yt + 1.5);

c = colorbar;
c.FontSize = 32;

numTicks = 6;
minVal = min(c.Limits);
maxVal = max(c.Limits);
ticks = linspace(minVal, maxVal, numTicks);

c.Ticks = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
c.FontSize = 32;

fprintf('MAGNITUDE_MEAN Min: %.3f\n', minVal);
fprintf('MAGNITUDE_MEAN Max: %.3f\n', maxVal);

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
    ax.XLim = [0 9];
    ax.XTickMode = 'manual';
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
    ax.YLim = [-1.5 1.5];
    ax.YTickMode = 'manual';
    newY = -1.5:1:1.5;
end
ax.YTick = newY;
ax.YTickLabel = arrayfun(@(y) sprintf('%d', y + 1.5), newY, 'UniformOutput', false);
% ===========================================================

saveas(fig_magMean, fullfile(contourf_folder, 'magnitude_mean_contour.png'));