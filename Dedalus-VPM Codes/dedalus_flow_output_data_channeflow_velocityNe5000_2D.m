clear all;
close all;
clc;

%folderpath = 'C:/Users/jinog/Documents/dedalus/channelflow/snapshots_channel'

% large channel box
%folderpath = 'E:/dedalus/Wavy Wall/dedalus_34241638_Re75/snapshots_channel'
folderpath = 'E:/dedalus/Wavy Wall/flexible signed distance function/opt_erf_peclet/dedalus_36868890__Re190_nek5000_hill_ubulk_2D/snapshots_channel'
%folderpath = 'E:/dedalus/Wavy Wall/flexible signed distance function/opt_erf_peclet/dedalus_local_Re100_nek5000_hill_2a_ubulk_2D/snapshots_channel'


%E:\dedalus\Wavy Wall\flexible signed distance function\opt_erf_peclet\dedalus_local_Re100_nek5000_hill_2c\snapshots_channel
%folderpath = 'C:/Users/jinog/Documents/dedalus/channelflow/checkpoints'
%C:\Users\jinog\Documents\dedalus\channelflow\checkpoints
cd(folderpath)

% dedalus_local_Re100_nek5000_hill_2a_ubulk_Re100
%folderpath = 'E:/dedalus/Wavy Wall/dedalus_33622799_Re75/snapshots_channel'
addpath('C:/Users/jinog/Documents/MATLAB/cbrewer2-master/cbrewer2');

%addpath('E:\dedalus\Wavy Wall\dedalus_32893242\snapshots_channel\cbrewer2-master');
%E:\dedalus\Wavy Wall\dedalus_local_Re75_3\snapshots_channel
% E:\dedalus\Wavy Wall\dedalus_34241638\snapshots_channel

% small channel box
%folderpath = 'E:/dedalus/small channel box/dedalus_15992631/snapshots_channel_stress'

%E:\dedalus\Large Channel Flow Box\dedalus_29471912

%folderpath = 'C:/Users/jinog/OneDrive - University of Connecticut/Documents/dedalus_29440202/snapshots_channel_stress'

% C:\Users\jinog\OneDrive - University of Connecticut\Documents\dedalus_29440202

% 29232568 _29277499
% C:\Users\jinog\Documents\dedalus\channelflow\dedalus_14148871 29376698

% C:\Users\jinog\OneDrive - University of Connecticut\Documents\dedalus_29372383\snapshots_channel_stress


% E:\dedalus\Wavy Wall\dedalus_34220049_Re75\snapshots_channel_stress
filePattern = fullfile(folderpath,'*.h5');
fileList = dir(filePattern);

flow = "2D" ; % flow can be 2D or 3D
%cd('C:/Users/jinog/Documents/MATLAB/snapshots') % cd into path where the file3s exist

h5disp('snapshots_channel_s2.h5'); % enter path to get to the file
%h5disp('checkpoints_s1.h5'); % enter path to get to the file

% folderpath2 = 'E:/dedalus/Wavy Wall/dedalus_34220049_Re75/snapshots_channel_stress'
% cd(folderpath2)
% h5disp('snapshots_channel_stress_s1.h5'); % enter path to get to the file


% Wavy validation case Choo parameters
epsilon = 0.5;
%epsilon = 0.25;
h = 1;
y0 = h;
% A1 = 0.54;
% A2 = 0.54;

A1 = epsilon;
A2 = epsilon;
%Lx = 0.6 * pi;
Lx = 9.0;
Ly = 1.0 + epsilon;

% For RE = 190
Ny = 576; % 48, 80
Nx = 384; % 36, 64

% For RE = 100
%Ny = 512; % 48, 80
%Nx = 352; % 36, 64
Nz = 2;

x_val = linspace(0, Lx, Nx);  % replaces fourdif x values for plotting only
y1_vals =  -y0 + A1*tanh(3.5*abs(x_val-4.5)- 3.5)

% y2_vals =  y0 + A2 * sin(2 * pi * x_val / Lx); % top wall
%mask_const = 75;
%% ------------------------------------------------------------------------
% physical properties
% nu = 1*10^-3 ; % dynamic viscosity of water
% h = 2; % height of channel
rho = 1000; % density of water

% Detect and decide Re based on file read
tok = regexp(folderpath, 'Re(\d+(_\d+)*)', 'tokens', 'once');
re_str = strrep(tok{1}, '_', '.');   % → '32.1'
% Convert to a number
Re = str2double(re_str);
% Now Re = 32.1 and you can use it throughout your script
fprintf('Detected Re = %.2f\n', Re);

show_velocity = "false";
mask_axis = "false"
video = "OFF"
data_t_tot = zeros(400,1);
n= 1;

%% --- discover the x_/y_ dataset names -----------------------
firstFile = fullfile(folderpath, fileList(1).name);
scale_info = h5info(firstFile, '/scales');
names      = {scale_info.Datasets.Name};

ix = find(startsWith(names,'x_'),1);
iy = find(startsWith(names,'y_'),1);
if isempty(ix)||isempty(iy)
    error('Couldn’t find x_* or y_* under /scales');
end

x_name = names{ix};
y_name = names{iy};
%% ------------------xxxxxxxxxxxxxxxxxx------------------------

% delete(gcp("nocreate"));
% parpool(8); % EQUAL TO --NTASKS IN CLUSTER

for k = 3:4 %length(fileList)

    %filePath = fullfile(folderpath, fileList(k).name);


    filePath = [folderpath,'/snapshots_channel_s',num2str(k),'.h5']

    % fprintf('Reading file: %s\n', filePath);
    %%%--------------------------------------------------------%%%%%
    dataV{k} = h5read(filePath, '/tasks/velocity');
    %%%--------------------------------------------------------%%%%%

    %dataK{k} = h5read(filePath, '/tasks/stiffness');

    %data1 = h5read(filePath, '/scales/sim_time');
    %data7{k} = h5read(filePath, '/tasks/stiffness');

    %% Velocity gradient
    dUdx{k} = h5read(filePath, '/tasks/dUdx');
    dUdy{k} = h5read(filePath, '/tasks/dUdy');
    dVdx{k} = h5read(filePath, '/tasks/dVdx');
    dVdy{k} = h5read(filePath, '/tasks/dVdy');

    %if k ==1
    %data_y = h5read(filePath, '/scales/y_hash_2c1da5940ac8593d53863fe7b647e57b44b9820f');
    %data_x = h5read(filePath, '/scales/x_hash_d0f963158a141cd8c7f416f18cfff03a6d595436');

    data_x = h5read(filePath, ['/scales/' x_name]);
    data_y = h5read(filePath, ['/scales/' y_name]);

    data_time{k} = h5read(filePath, '/scales/sim_time');
    data_mask{k} = h5read(filePath, '/tasks/mask');

    %end
    k;

end

%% Call function to create mask
main_signed_dist_functn_nek5000_1a;

% y1_vals = data_y;
% x_val = data_x;
size(dataV)
%load(dataV)

dataV_all = cat(4, dataV{:});
dUdx_all = cat(3, dUdx{:});
dUdy_all = cat(3, dUdy{:});
dVdx_all = cat(3, dVdx{:});
dVdy_all = cat(3, dVdy{:});

v1 = squeeze(dataV_all(:,:,1,:));
v2 = squeeze(dataV_all(:,:,2,:));

% U_bar = mean(U_bar(:,end-4000:end),2);

v1_avg = mean(v1,3);
v2_avg = mean(v2,3);
u_mean_zt = v1_avg;
v_mean_zt = v2_avg;

v1_avg(mask_smooth_solid) = NaN;
v2_avg(mask_smooth_solid) = NaN;
magnitude_mean = sqrt(v1_avg.^2 + v2_avg.^2)

dUdx_mean_zt = mean(dUdx_all,3);
dUdy_mean_zt = mean(dUdy_all,3);
dVdx_mean_zt = mean(dVdx_all,3);
dVdy_mean_zt = mean(dVdy_all,3);
mask_smooth = squeeze(data_mask{k}(:,:,100));
% mask_smooth = permute(mask_smooth,[2 1 3]);

% v1_avg = permute(v1_avg,[3,2,1]);
% v2_avg = permute(v2_avg,[3,2,1]);
% v3_avg = permute(v3_avg,[3,2,1]);
%
% dUdx_avg = permute(dUdx_avg,[3,2,1]);
% dUdy_avg = permute(dUdy_avg,[3,2,1]);
% dVdx_avg = permute(dVdx_avg,[3,2,1]);
% dVdy_avg = permute(dVdy_avg,[3,2,1]);

% make output folder with Re in the name, e.g. mean_Re75.30
outdir = fullfile(pwd, sprintf('mean_v_Re%.2f_c%.2f', Re,c));
if ~exist(outdir,'dir'), mkdir(outdir); end

% save each to its own .mat
save(fullfile(outdir,'u_mean_zt.mat'),'u_mean_zt');       % add '-v7.3' if >2GB
save(fullfile(outdir,'v_mean_zt.mat'),'v_mean_zt');


save(fullfile(outdir,'dUdx_mean_zt.mat'),'dUdx_mean_zt');       % add '-v7.3' if >2GB
save(fullfile(outdir,'dUdy_mean_zt.mat'),'dUdy_mean_zt');
save(fullfile(outdir,'dVdx_mean_zt.mat'),'dVdx_mean_zt');
save(fullfile(outdir,'dVdy_mean_zt.mat'),'dVdy_mean_zt');

save(fullfile(outdir,'data_x.mat'),'data_x');
save(fullfile(outdir,'data_y.mat'),'data_y');
save(fullfile(outdir,'mask_smooth.mat'),'mask_smooth');
%% ------------------------------------------------------------------------
%% -----------------------------END ---------------------------------------

ntime = size(dataV_all, 4);   % total number of time steps
Nx = size(dataV_all, 2)
Ny = size(dataV_all, 1)
nsnaps = 30;               % number of snapshots you want

% Time indices to extract (evenly spaced)
snap_idx = round(linspace(2, ntime, nsnaps));

% Preallocate cell arrays
u_xy = cell(1, nsnaps);
v_xy = cell(1, nsnaps);
u_yz = cell(1, nsnaps);
v_yz = cell(1, nsnaps);

for i = 1:nsnaps
    % Extract and permute velocity components
    u = squeeze(dataV_all(:,:,1,snap_idx(i)));
    v = squeeze(dataV_all(:,:,2,snap_idx(i)));


    % Extract xy and yz slices
    u_xy{i} = squeeze(u(:,:));
    v_xy{i} = squeeze(v(:,:));

end

[X, Y] = meshgrid(data_x, data_y);


nLevels = 30;
nLevelsv = 30;
tick_enforcing = "manual";
contourf_folder = fullfile("figures_snapshots", "contourf");
if show_velocity == "true"

    % Create subfolders
    %pcolor_folder = fullfile("figures_snapshots", "pcolor");
    %contour_folder = fullfile("figures_snapshots", "contour");
    contourf_folder = fullfile("figures_snapshots", "contourf");

    if ~exist(pcolor_folder, 'dir')
        mkdir(pcolor_folder);
    end
    if ~exist(contour_folder, 'dir')
        mkdir(contour_folder);
    end
    if ~exist(contourf_folder, 'dir')
        mkdir(contourf_folder);
    end
    rwb1 = bluewhitered;
    %% === PCOLOR FIGURE LOOP ===
    % for i = 1:nsnaps
    %     % --- u_xy pcolor ---
    %     fig1 = figure(2*(i-1) + 1);
    %     set(fig1, 'Position', [100, 100, 1000, 800]);
    %     pcolor(X, Y, u_xy{i});
    %     shading interp
    %     clim([min(u_xy{i}(:)) max(u_xy{i}(:))])
    %     colormap(rwb1);
    %
    %     c = colorbar;
    %     existingTicks = c.Ticks;
    %     minVal = min(c.Limits);
    %     maxVal = max(c.Limits);
    %     % Add min and max if not already present
    %     numTicks = 6;
    %     newTicks = unique([existingTicks, minVal, maxVal]);
    %     ticks = linspace(minVal, maxVal, numTicks);
    %     c.Ticks = ticks; % or newticks if you want max value
    %     % Update labels to show actual values
    %     c.TickLabels = arrayfun(@(x) num2str(x), ticks, 'UniformOutput', false);
    %     c.FontSize = 32;
    %
    %     %title(sprintf("U Velocity (pcolor) - Snapshot %d (t = %d)", i, snap_idx(i)))
    %     hold on
    %     %plot(x_val, y1_vals, 'k', 'LineWidth', 1.5)
    %     hold on
    %     plot(x_solid,y_solid,'.');
    %     hold on
    %     set(gca, 'FontSize', 32);
    %     xlabel('X', 'FontSize', 24);
    %     ylabel('Y', 'FontSize', 24);
    %     saveas(fig1, fullfile(pcolor_folder, sprintf('u_pcolor_snapshot_%d.png', i)));
    %
    %     % --- v_xy pcolor ---
    %     fig2 = figure(2*(i-1) + 2);
    %     set(fig2, 'Position', [100, 100, 1000, 800]);
    %     pcolor(X, Y, v_xy{i});
    %     shading interp
    %     clim([min(v_xy{i}(:)) max(v_xy{i}(:))])
    %     colormap(rwb1);
    %
    %     c = colorbar;
    %     existingTicks = c.Ticks;
    %     minVal = min(c.Limits);
    %     maxVal = max(c.Limits);
    %     % Add min and max if not already present
    %     newTicks = unique([existingTicks, minVal, maxVal]);
    %     ticks = linspace(minVal, maxVal, numTicks);
    %     c.Ticks = ticks;
    %     % Update labels to show actual values
    %     c.TickLabels = arrayfun(@(x) num2str(x), ticks, 'UniformOutput', false);
    %     c.FontSize = 32;
    %
    %     %title(sprintf("V Velocity (pcolor) - Snapshot %d (t = %d)", i, snap_idx(i)))
    %     hold on
    %     %plot(x_val, y1_vals, 'k', 'LineWidth', 1.5)
    %     hold on
    %     plot(x_solid,y_solid,'.');
    %     hold on
    %     set(gca, 'FontSize', 32);
    %     xlabel('X', 'FontSize', 24);
    %     ylabel('Y', 'FontSize', 24);
    %     saveas(fig2, fullfile(pcolor_folder, sprintf('v_pcolor_snapshot_%d.png', i)));
    % end
    %
    % clf(fig1)
    % clf(fig2)
    % %% === CONTOUR FIGURE LOOP ===
    % for i = 1:nsnaps
    %     % --- u_xy contour ---
    %     fig3 = figure(2*(i-1) + 1);
    %     clf(fig3)
    %     set(fig3, 'Position', [100, 100, 1000, 800]);
    %     contour(X, Y, u_xy{i}, nLevels + 20);  % Removed LineStyle
    %     colormap(rwb1);
    %
    %       % === Y-AXIS RELABELING (+1.5 shift) ===
    %     yt = get(gca, 'YTick');                    % original ticks (e.g., -1.5..1.5)
    %     set(gca, 'YTickLabel', yt + 1.5);          % relabeled as 0..3
    %
    %     c = colorbar;
    %     existingTicks = c.Ticks;
    %     minVal = min(c.Limits);
    %     maxVal = max(c.Limits);
    %     numTicks = 6;
    %     % Add min and max if not already present
    %     newTicks = unique([existingTicks, minVal, maxVal]);
    %     ticks = linspace(minVal, maxVal, numTicks);
    %     c.Ticks = ticks;
    %     % Update labels to show actual values
    %     c.TickLabels = arrayfun(@(x) num2str(x), ticks, 'UniformOutput', false);
    %     c.FontSize = 24;
    %
    %     title(sprintf("U Velocity (contour) - Snapshot %d (t = %d)", i, snap_idx(i)))
    %     hold on
    %     %plot(x_val, y1_vals, 'k', 'LineWidth', 1.5)
    % %     hold on
    %     plot(x_solid,y_solid,'.');
    %     hold off
    %     set(gca, 'FontSize', 32);
    %     xlabel('X', 'FontSize', 24);
    %     ylabel('Y', 'FontSize', 24);
    %     saveas(fig3, fullfile(contour_folder, sprintf('u_contour_snapshot_%d.png', i)));
    %
    %     % --- v_xy contour ---
    %     fig4 = figure(2*(i-1) + 2);
    %     clf(fig4)
    %     set(fig4, 'Position', [100, 100, 1000, 800]);
    %     contour(X, Y, v_xy{i}, nLevelsv + 20);  % Removed LineStyle
    %     colormap(rwb1);
    %
    %       % === Y-AXIS RELABELING (+1.5 shift) ===
    %     yt = get(gca, 'YTick');                    % original ticks (e.g., -1.5..1.5)
    %     set(gca, 'YTickLabel', yt + 1.5);          % relabeled as 0..3
    %
    %     c = colorbar;
    %     existingTicks = c.Ticks;
    %     minVal = min(c.Limits);
    %     maxVal = max(c.Limits);
    %     % Add min and max if not already present
    %     newTicks = unique([existingTicks, minVal, maxVal]);
    %     ticks = linspace(minVal, maxVal, numTicks);
    %     c.Ticks = ticks;
    %     % Update labels to show actual values
    %     c.TickLabels = arrayfun(@(x) num2str(x), ticks, 'UniformOutput', false);
    %     c.FontSize = 24;
    %
    %     title(sprintf("V Velocity (contour) - Snapshot %d (t = %d)", i, snap_idx(i)))
    %     hold on
    %     %plot(x_val, y1_vals, 'k', 'LineWidth', 1.5)
    % %     hold on
    %     plot(x_solid,y_solid,'.');
    %     hold off
    %     set(gca, 'FontSize', 32);
    %     xlabel('X', 'FontSize', 24);
    %     ylabel('Y', 'FontSize', 24);
    %     saveas(fig4, fullfile(contour_folder, sprintf('v_contour_snapshot_%d.png', i)));
    % end

    tick_enforcing = "manual";
    %% === CONTOURF FIGURE LOOP ===
    for i = 1:nsnaps
        % --- u_xy contour ---
        fig5 = figure(2*(i-1) + 1);
        clf(fig5)
        set(fig5, 'Position', [100, 100, 1000, 800]);
        contourf(X, Y, u_xy{i}, nLevels);
        colormap(rwb1);
        %     xticks(linspace(0,Lx,9))
        %     yticks(linspace(data_y(1),data_y(end),6))

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

        saveas(fig5, fullfile(contourf_folder, sprintf('u_contour_snapshot_%d.png', i)));

        %% --- v_xy contour --- %%
        fig6 = figure(2*(i-1) + 2);
        clf(fig6)
        set(fig6, 'Position', [100, 100, 1000, 800]);
        contourf(X, Y, v_xy{i}, nLevelsv);
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
        saveas(fig6, fullfile(contourf_folder, sprintf('v_contour_snapshot_%d.png', i)));

    end

    %mean_velocity_dedalus_postprocess;


end
rwb1 = bluewhitered;
mean_velocity_dedalus_postprocess;

if video == "ON"
    create_video_nekhill;
end
u_dedalus = u_xy{12};
v_dedalus = v_xy{12};

save('u_dedalus.mat','u_xy')
save('v_dedalus.mat','v_xy')

save('dedalus_x.mat','data_x')
save('dedalus_y.mat','data_y')



%% ------------------------xxxxxxxxxxxxxxxxxxxxxxxxxx----------------------
%% Find u velocity slices at given x to understand no-slip location
%% ------------------------xxxxxxxxxxxxxxxxxxxxxxxxxx----------------------
% Create output directory for plots
output_dir = sprintf('plots_Re_%g_mask_%.1f', Re, mask_const);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Pick centerline slices
n = Nx;
x_idx_vec = linspace(1,n,n);  % <-- EDIT THIS VECTOR
z_idx = 1;          % center in z
t_idx = snap_idx(end);        % last snapshot
u_thresh = 1e-3;              % velocity threshold
x_tol = 1e-6;                 % tolerance for wall match
%mask_smooth = squeeze(mask_smooth(:,:,t_idx));

v1 = permute(v1, [2, 1, 3]);  % reorder to (x, y, z, t)
v2 = permute(v2, [2, 1, 3]);  % reorder to (x, y, z, t)
%mask_smooth = mask_smooth';     % transpose if needed
l_star_grid = nan(length(x_idx_vec), 2);  % preallocate (2 l* per x)

%mask_midline_hillp;

%ywall_values;
mask_smooth = permute(mask_smooth,[2,1]);

for k = 1:length(x_idx_vec)
    x_idx = x_idx_vec(k);

    % Extract velocity and mask
    u_y = squeeze(v1(x_idx, :, t_idx));
    v_y = squeeze(v2(x_idx, :, t_idx));
    mask_y = squeeze(mask_smooth(x_idx, :));

    %mask_y_calc;

    % Match x to wavy wall location
    %x_match = abs(x_wall2 - data_x(x_idx)) < x_tol;
    x_match = data_x(x_idx);
    y_wall_here1 = y_wall2a(x_idx);
    %y_wall_here2 = y_wall2b(x_idx);
    y_wall_here = [y_wall_here1];

    [~, iy1] = min(abs(data_y - y_wall_here(1)));
    %  [~, iy2] = min(abs(data_y - y_wall_here(2)));

    y_wall_here1_cheb = data_y(iy1);   % scalar (closest y-grid point to bottom wall)
    % y_wall_here2_cheb = data_y(iy2);   % scalar (closest y-grid point to top wall)

    if mod(x_idx,40) == 0
        %% ----------------Plot velocity(u,v) and mask ------------------------
        %% --------------------u velocity vs mask------------------------------
        figure;
        yyaxis left
        plot(data_y, u_y, 'b', 'LineWidth', 2);
        ylabel('Velocity U','FontSize', 20);
        %ylim([min(u_y)- 0.25, max(u_y)])
        ylim([min(u_y)- 0.25, 3.0])

        if mask_axis == "true"
            yyaxis right
            plot(data_y, mask_y, 'ko', 'LineStyle','none', 'MarkerFaceColor','k', 'MarkerSize',4);
            ylabel('Mask','FontSize', 20);
            %ylim([0, max(mask_y)*1.1])
            ylim([0 120]);   % <<< FIXED MASK RANGE
        end
        %xlabel('y (wall normal)','FontSize', 20);
        xlabel('y','FontSize', 20);
        %title(sprintf('U Velocity Profile vs. Mask at x = %.3f, t = %d, mask = %.1f', data_x(x_idx), t_idx, mask_const));
        grid on;


        ytick_range = min(data_y):1.0:3.0;
        xticks(ytick_range);
        xticklabels(arrayfun(@(x) sprintf('%.1f', x+1.5), ytick_range, 'UniformOutput', false));
        %xlim([ytick_range(1) ytick_range(end)]);

        % --- MAKE TICKS LARGER (THE ONLY CHANGE YOU WANTED) ---
        ax = gca;
        ax.XAxis.FontSize = 24;          % x ticks
        ax.YAxis(1).FontSize = 24;       % left y ticks
        ax.YAxis(2).TickValues = [];    % remove ticks
        ax.YAxis(2).TickLabels = {};    % remove tick labels

        % ==========================================================
        %         ADD MASK MIDLINE OVERLAY (from mask_y)
        % ==========================================================
        mask_norm = mask_y / max(mask_y);
        sgn = mask_norm - 0.5;
        idx = find(sgn(1:end-1) .* sgn(2:end) <= 0, 1, 'first');

        if ~isempty(idx)
            y1 = data_y(idx);      y2 = data_y(idx+1);
            m1 = mask_norm(idx);   m2 = mask_norm(idx+1);
            mask_y_midline = y1 + (0.5 - m1) * (y2 - y1) / (m2 - m1);

            yyaxis left
            xline(mask_y_midline, 'r--', 'LineWidth', 2);
            % optional label:

        end

        if mask_axis == "true"
            lgd = legend('U(x,y)','y_{hill}', 'mask(y)', 'Location', 'northeast');
        else
            lgd = legend('U(x,y)','y_{hill}', 'Location', 'northeast');
        end
        lgd.FontSize = 24;
        lgd.Color = 'none';   % makes background transparent
        lgd.Box = 'off';
        % ==========================================================
        filename_safe = sprintf('u velocity_mask_x_%.3f_t_%d_mask_%.1f.png', data_x(x_idx), t_idx, mask_const);
        saveas(gcf, fullfile(output_dir, filename_safe));


        %% end of plot --------------xxxxxxxxxxxxxxxxxx------------------------
    end





    %% --------------Displacement length calculation-----------------------
    %     crossing_idx_up   = find((u_y(1:end-1) < u_thresh) & (u_y(2:end) > u_thresh));
    %     crossing_idx_down = find((u_y(1:end-1) > u_thresh) & (u_y(2:end) < u_thresh));
    %     crossing_idx = sort([crossing_idx_up; crossing_idx_down]);


    %     crossing_idx_up = [];
    %     crossing_idx_down = [];
    %
    %     % Search for crossings with early exit
    %     for i = 1:length(u_y)-1
    %         if isempty(crossing_idx_up) && u_y(i) <= u_thresh && u_y(i+1) >= u_thresh
    %             crossing_idx_up = i;
    %         end
    %         if isempty(crossing_idx_down) && u_y(i) >= u_thresh && u_y(i+1) <= u_thresh
    %             crossing_idx_down = i;
    %         end
    %         if ~isempty(crossing_idx_up) && ~isempty(crossing_idx_down)
    %             break;
    %         end
    %     end
    %
    %     % Combine in order (bottom wall first)
    %     %crossing_idx = sort([crossing_idx_up, crossing_idx_down]);
    %     crossing_idx = sort(crossing_idx_up);
    %
    %     l_star = zeros(1, length(crossing_idx));
    %     y_start = zeros(1, length(crossing_idx));
    %
    %     for j = 1:length(crossing_idx)
    %         i = crossing_idx(j);
    %         y1 = data_y(i); y2 = data_y(i+1);
    %         u1 = u_y(i);    u2 = u_y(i+1);
    %         y_start(j) = y1 + (u_thresh - u1) * (y2 - y1) / (u2 - u1);
    %         l_star(j) = y_start(j) - y_wall_here(j);
    %     end
    %
    %
    %     if length(l_star) == 1
    %         l_star_grid(k) = l_star;        % scalar per x-index
    %     else
    %         warning('At x = %.3f: expected 1 crossing, got %d', ...
    %             data_x(x_idx), length(l_star));
    %         l_star_grid(k) = nan;
    %     end
    %
    %
    %     l_star_bottom = l_star(1);
    %     y_wall_bottom = y_wall_here(1);
    %

end


%
%     T = table(l_star_grid(:,1), l_star_grid(:,2), ...
%     'VariableNames', {'l_star_lower', 'l_star_upper'});
%     writetable(T, fullfile(output_dir, 'l_star_grid.csv'));