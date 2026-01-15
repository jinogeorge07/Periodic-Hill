% Read Nek5000 output data
clc; clear all;close all;

folderpath = './'

hillp_out = readmatrix("hilllp_refined_Re100.csv");

%hillp_output_Re100_timestep_1840
%C:\Users\jinog\Documents\MATLAB\hillp

addpath('C:/Users/jinog/Documents/MATLAB/cbrewer2-master/cbrewer2');

nx = 64;
ny = 309;
nlevels = 30;

% columns in hillp_out:
% 1=u, 2=v, ..., 6=x, 7=y
x_all = hillp_out(:,6);
y_all = hillp_out(:,7);
u_all = hillp_out(:,1);
v_all = hillp_out(:,2);

x_list = unique(x_all);

% how many y-samples per x-column (some columns may be oversampled)
counts = zeros(numel(x_list),1);
for i = 1:numel(x_list)
    counts(i) = sum(x_all == x_list(i));
end
base_count = mode(counts);   % baseline Ny per column (e.g., 64)

% (optional) keep your ny consistent with what's in the file
% if you already set ny, you can assert they match or just use base_count.
% ny = base_count;

% preallocate using the detected baseline
y_2D = nan(numel(x_list), base_count);
x_2D = nan(numel(x_list), base_count);
u_2D = nan(numel(x_list), base_count);
v_2D = nan(numel(x_list), base_count);

for xi = 1:numel(x_list)
    nek_ind = find(x_all == x_list(xi));
    % sort this column by y
    col = sortrows([y_all(nek_ind), u_all(nek_ind), v_all(nek_ind), x_all(nek_ind)], 1);
    ny_here = size(col,1);

    if ny_here == base_count
        take = 1:ny_here;
    else
        % integer decimation if oversampled (2x, 3x, ...)
        k = max(1, round(ny_here / base_count));
        take = 1:k:ny_here;
        % enforce exact length (handles tiny non-integer ratios / dup nodes)
        if numel(take) ~= base_count
            take = round(linspace(1, ny_here, base_count));
        end
    end

    y_2D(xi,:) = col(take,1).';
    u_2D(xi,:) = col(take,2).';
    v_2D(xi,:) = col(take,3).';
    x_2D(xi,:) = col(take,4).';
end

rwb1 = bluewhitered;
%% U nek 
fig = figure;
%set(fig, 'Position', [100, 100, 900, 750]);  % Resize figure
set(fig, 'Position', [100, 100, 1000, 800]); 
contourf(x_2D, y_2D, u_2D, nlevels);
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
saveas(gcf, 'U_nek.png');


%% V nek ------------------------------------------------------------------
fig = figure;
%set(fig, 'Position', [100, 100, 900, 750]);  % Resize figure
set(fig, 'Position', [100, 100, 1000, 800]); 
contourf(x_2D, y_2D, v_2D, nlevels);
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
saveas(gcf, 'V_nek.png');