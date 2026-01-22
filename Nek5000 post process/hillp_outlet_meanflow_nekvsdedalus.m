% Read Nek5000 output data
clc; clear all;close all;
folderpath = 'C:/Users/jinog/Documents/MATLAB/hillp'
cd(folderpath)


%hillp_out{1} = readmatrix("hilllp_refined_Re190_15.csv");

%hillp_output_Re100_timestep_1840
%C:\Users\jinog\Documents\MATLAB\hillp

addpath('C:/Users/jinog/Documents/MATLAB/cbrewer2-master/cbrewer2');

nx = 64;
ny = 309;
nlevels = 30;


% --- settings ---
file_start = 700;
file_end   = 700;
nf = file_end - file_start + 1;
rwb1 = bluewhitered;
% --- read first file as reference ---
for i=1:nf
    %fname{i} = ['nek_mean_flow_750-1500timeunits_U.csv'];
    fname{i} = ['nek_mean_flow_1500-2500timeunits.csv'];
    hillp_out{i} = readmatrix(fname{i});


    % columns in hillp_out:
    % 1=u, 2=v, ..., 6=x, 7=y
    % x_all = hillp_out{i}(:,5);
    % y_all = hillp_out{i}(:,6);
    % u_all = hillp_out{i}(:,1);
    % v_all = hillp_out{i}(:,2);


    x_all = hillp_out{i}(:,4);
    y_all = hillp_out{i}(:,5);
    u_all = hillp_out{i}(:,1);
    v_all = hillp_out{i}(:,2);
    x_list = unique(x_all);

    % how many y-samples per x-column (some columns may be oversampled)
    counts = zeros(numel(x_list),1);
    for ii = 1:numel(x_list)
        counts(ii) = sum(x_all == x_list(ii));
    end
    base_count = mode(counts);   % baseline Ny per column (e.g., 64)

    % (optional) keep your ny consistent with what's in the file
    % if you already set ny, you can assert they match or just use base_count.
    % ny = base_count;

    % preallocate using the detected baseline
    y_2D{i} = nan(numel(x_list), base_count);
    x_2D{i} = nan(numel(x_list), base_count);
    u_2D{i} = nan(numel(x_list), base_count);
    v_2D{i} = nan(numel(x_list), base_count);

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

        y_2D{i}(xi,:) = col(take,1).';
        u_2D{i}(xi,:) = col(take,2).';
        v_2D{i}(xi,:) = col(take,3).';
        x_2D{i}(xi,:) = col(take,4).';
    end

    % figure;
    % contourf(x_2D{i},y_2D{i},u_2D{i},nlevels);
    % shading interp
    % colormap(rwb1);
    % colorbar
    % hold on
    %
    % figure;
    % contourf(x_2D{i},y_2D{i},v_2D{i},nlevels);
    % shading interp
    % colormap(rwb1);
    % colorbar
    % hold on

end

% ---- STACK + MEAN (assumes all grids have same size/order) ----
Nx = size(u_2D{1},1);
Ny = size(u_2D{1},2);

Ustack = nan(Nx,Ny,nf);
Vstack = nan(Nx,Ny,nf);

for i = 1:nf
    %     if ~isequal(size(u_2D{i}), [Nx,Ny]) || ~isequal(size(v_2D{i}), [Nx,Ny])
    %         error('Size mismatch at i=%d. Some timestep has different base_count/Nx.', i);
    %     end
    %     % Optional but recommended: ensure grids match too
    %     if max(abs(x_2D{i}(:) - x_2D{1}(:))) > 1e-12 || max(abs(y_2D{i}(:) - y_2D{1}(:))) > 1e-12
    %         error('Grid mismatch at i=%d (x_2D/y_2D differs).', i);
    %     end

    u_2D_mean = u_2D{i};
    v_2D_mean = v_2D{i};
end


% Use the grid from the first file
x_mean = x_2D{1};
y_mean = y_2D{1};

%% U mean nek

fig = figure;
set(fig, 'Position', [100, 100, 1000, 800]);
contourf(x_mean, y_mean, u_2D_mean, nlevels);
shading interp
colormap(rwb1);

% -------- Colorbar setup (unchanged, but cleaned) ----------
c = colorbar;
c.Limits = [min(u_2D_mean(:)) max(u_2D_mean(:))];
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

saveas(gcf, 'U_nek.png');



%% V nek ------------------------------------------------------------------
fig = figure;
set(fig, 'Position', [100, 100, 1000, 800]);

contourf(x_mean, y_mean, v_2D_mean, nlevels);
shading interp
colormap(rwb1);

% ---------------- COLORBAR ----------------
c = colorbar;
c.Limits = [min(v_2D_mean(:)) max(v_2D_mean(:))];
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

saveas(gcf, 'V_nek.png');
%% ------------------------------------------------------------------------

x_2D_new = x_2D{end};
y_2D_new = y_2D{end};


x_flat = x_2D_new(:);
y_flat = y_2D_new(:);
u_flat = u_2D_mean(:);
v_flat = v_2D_mean(:);

y_flat = y_flat - 1.5;



%% Dedalus output
%load('u_dedalus.mat')
%load('v_dedalus.mat')

load('u_mean_zt.mat')
load('v_mean_zt.mat')
load('data_x.mat')
load('data_y.mat')
load('X_fluid.mat')
load('Y_fluid.mat')

[X,Y] = meshgrid(data_x,data_y);


u_xy1 = u_xy{12};
v_xy1 = v_xy{12};

% % Evaluate on Dedalus grid
% F_u = griddata(x_flat,y_flat,u_flat,x_new,y_new);
% F_v = griddata(x_flat,y_flat,v_flat,x_new,y_new);
%
% F_u = reshape(F_u,size(X));
% F_v = reshape(F_v,size(Y));

%U_nek2ded = interp2(x_2D,y_2D,u_2D,X_fluid,Y_fluid);

% Create interpolants
F_u = scatteredInterpolant(x_flat, y_flat, u_flat, 'linear',  'nearest');
F_v = scatteredInterpolant(x_flat, y_flat, v_flat, 'linear',  'nearest');

U_nek2ded = F_u(X_fluid, Y_fluid);

V_nek2ded = F_v(X_fluid, Y_fluid);

% U_nek2ded_ts1840 = F_u(X,Y);
% V_nek2ded_ts1840 = F_v(X,Y);

% save('u_nek2ded_ts2000.mat','U_nek2ded')
% save('v_nek2ded_ts2000.mat','V_nek2ded')

% load("u_nek2ded_ts2000.mat")
% load("v_nek2ded_ts2000.mat")
%
% diff_u_nek = U_nek2ded_ts1840 - U_nek2ded;
% diff_v_nek = V_nek2ded_ts1840 - V_nek2ded;
%
% % see if velocity has reached steady state
% figure;
% plot(diff_u_nek);
% title("u velocity has reached steady state")
% figure;
% plot(diff_v_nek);
% title("v velocity has reached steady state")

% figure;
% contour(X,Y,U_nek2ded,nlevels);
% shading interp
% colormap(rwb1);
% title(sprintf("U Velocity  interpolation from nek5000 to dedalus grid"))
% c = colorbar;
% existingTicks = c.Ticks;
% minVal = min(c.Limits);
% maxVal = max(c.Limits);
% % Add min and max if not already present
% newTicks = unique([existingTicks, minVal, maxVal]);
% c.Ticks = newTicks;
% % Update labels to show actual values
% c.TickLabels = arrayfun(@(x) num2str(x), newTicks, 'UniformOutput', false);
% saveas(gcf, 'U_nek2ded_interp1.png');
% hold on
%
%
% figure;
% contour(X,Y,V_nek2ded,nlevels);
% shading interp
% colormap(rwb1);
% title(sprintf("V Velocity  interpolation from nek5000 to dedalus grid"))
% c = colorbar;
% existingTicks = c.Ticks;
% minVal = min(c.Limits);
% maxVal = max(c.Limits);
% % Add min and max if not already present
% newTicks = unique([existingTicks, minVal, maxVal]);
% c.Ticks = newTicks;
% % Update labels to show actual values
% c.TickLabels = arrayfun(@(x) num2str(x), newTicks, 'UniformOutput', false);
% saveas(gcf, 'V_nek2ded_interp1.png');
% colorbar
% hold on

%% U nek to Dedalus ( curvilinear to cartesiasn grid)
fig = figure;
%set(fig, 'Position', [100, 100, 900, 750]);  % Resize figure
set(fig, 'Position', [100, 100, 1000, 800]);
contourf(X, Y, U_nek2ded, nlevels);
shading interp
colormap(rwb1);

% title("U Velocity interpolation from nek5000 to dedalus grid", 'FontSize', 22);

% Create and format colorbar
c = colorbar;
c.FontSize = 32;  % Increase colorbar tick label size

% Remove min/max from ticks
existingTicks = c.Ticks;
numTicks = 6;
c.Limits = [-0.15 1.93];
minVal = min(c.Limits);
maxVal = max(c.Limits);
ticks = linspace(minVal, maxVal, numTicks);% Remove min and max from ticks

% Filter out min/max from tick list
filteredTicks = existingTicks(existingTicks ~= minVal & existingTicks ~= maxVal);
% c.Ticks = filteredTicks;
c.Ticks = linspace(c.Limits(1), c.Limits(2), numTicks);
c.Ticks = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);

% Print min/max separately
fprintf('U velocity Nek to dedalus grid min Value: %.3f\n', minVal);
fprintf('U velocity Nek to dedalus grid xax Value: %.3f\n', maxVal);

% Axis formatting
set(gca, 'FontSize', 32);  % Axis tick label size
xlabel('X', 'FontSize', 24);
ylabel('Y', 'FontSize', 24);

saveas(gcf, 'U_nek2ded_interp2.png');
hold on

%% V nek to Dedalus ( curvilinear to cartesiasn grid)
fig = figure;
%set(fig, 'Position', [100, 100, 900, 750]);  % Resize figure
set(fig, 'Position', [100, 100, 1000, 800]);
contourf(X, Y, V_nek2ded, nlevels);
shading interp
colormap(rwb1);

% title("V Velocity interpolation from nek5000 to dedalus grid", 'FontSize', 22);

% Create and format colorbar
c = colorbar;
c.FontSize = 32;  % Increase colorbar tick label size
% Remove min/max from ticks
numTicks = 6;
existingTicks = c.Ticks;
c.Limits = [-0.07 0.28]
minVal = min(c.Limits);
maxVal = max(c.Limits);
ticks = linspace(minVal, maxVal, numTicks);% Remove min and max from ticks
filteredTicks = existingTicks(existingTicks ~= minVal & existingTicks ~= maxVal);
% c.Ticks = filteredTicks;
% c.Ticks = linspace(c.Limits(1), c.Limits(2), numTicks);
c.Ticks = ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);

% Print min/max separately
fprintf('V velocity Nek to dedalus grid min Value: %.3f\n', minVal);
fprintf('V velocity Nek to dedalus grid max Value: %.3f\n', maxVal);

% Axis formatting
set(gca, 'FontSize', 32);  % Axis tick label size
xlabel('X', 'FontSize', 24);
ylabel('Y', 'FontSize', 24);

saveas(gcf, 'V_nek2ded_interp2.png');
hold on
%% Error Computation

% Assumes you already have:
% U_nek2ded, V_nek2ded   % nek velocities mapped onto Dedalus grid
% u_xy1, v_xy1           % Dedalus velocities on same grid (note the transpose below)
% data_x (1×nx), data_y (ny×1)  % optional, for locating the max
% [X,Y] = meshgrid(data_x, data_y);  % optional

% 1) Valid mask using YOUR terms/orientation
mask = isfinite(U_nek2ded) & isfinite(V_nek2ded) & ...
    isfinite(u_xy1)    & isfinite(v_xy1);

% 2) Component errors (your orientation: Dedalus transposed)
error_x = U_nek2ded - u_xy1;
error_y = V_nek2ded - v_xy1;

% 3) Vector energy magnitude
E = (error_x.^2 + error_y.^2);

% 4) Basic stats on valid points
Error = E(mask);
err_min  = min(Error);
err_max  = max(Error);
err_mean = mean(Error);
err_rmse = sqrt(mean(Error.^2));
err_mae  = mean(abs(Error));

fprintf('Vector error |Δu|:  min=%.4e  mean=%.4e  RMSE=%.4e  max=%.4e\n', ...
    err_min, err_mean, err_rmse, err_max);

% 5) Relative errors (vs Nek magnitude)

Sref = sqrt(mean(mean(U_nek2ded(mask).^2 + V_nek2ded(mask).^2)));
Rel_rmse  = err_rmse ./ Sref;

rel_u = error_x/max(u_xy1(:));
rel_v = error_y/max(v_xy1(:));

fprintf('Relative to |u_nek|:  RMSE=%.4e\n', Rel_rmse);

% 6) Where is the max error?
E_masked = E; E_masked(~mask) = -Inf;           % exclude invalids
[Em, linIdx] = max(E_masked(:));
[iMax, jMax] = ind2sub(size(E_masked), linIdx);
fprintf('Max at indices (i=%d, j=%d)\n', iMax, jMax);

% If you have X,Y from meshgrid(data_x,data_y):
try
    fprintf('Max at (x=%.6f, y=%.6f)\n', X(iMax,jMax), Y(iMax,jMax));
catch
    % omit if X,Y not defined
end

% 7) Quick visualization (optional)
figure; contourf(E,30,'LineStyle','none'); colorbar
set(gca,'YDir','normal'); axis tight equal
title('|u_{ded} - u_{nek}| (vector error magnitude)');
saveas(gcf, 'error.png');
hold on; plot(jMax, iMax, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5); hold off