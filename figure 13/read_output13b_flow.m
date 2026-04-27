% Input Output Post Processing for 2 way coupled flow
clc; clear all; close all;
%% Add path for color profile---------------------------------------------------------
%=========================================================================================================================
folderpath = './'

%% Periodic hill and Domain Properties
Ny = 168; Nx = 120; Nz = 2; % Nz to Nx for streamwise terms
N = Ny*Nx;
epsilon = 0.50;
%epsilon = 0.25;
h = 1;
y0 = h;
A1 = epsilon;
A2 = epsilon;
Lx = 9.0;
Ly = 1.0 + 1.0*epsilon;
Re = 300
x_val = linspace(0, Lx, Nx);  % replaces fourdif x values for plotting only
A  = 4.5;
B  = 3.5;
C  = 1/6;

% Read the 3 files directly and stack the already-saved fields
S1 = load('stability_results_Re300_168x120_kz20_c1.mat');
S4 = load('stability_results_Re300_168x120_kz20_c4.mat');
S6 = load('stability_results_Re300_168x120_kz20_c6.mat');

% take x, y, kz info from one file
x  = S1.x;
y  = S1.y;
kz_idx = 20;
c_index = [1 4 6];
nc = numel(c_index);
kz = 0.55;
% pick only the desired c-slice from each file
U_hat  = cat(3, S1.U_hat(:,:,1),  S4.U_hat(:,:,4),  S6.U_hat(:,:,6));
V_hat  = cat(3, S1.V_hat(:,:,1),  S4.V_hat(:,:,4),  S6.V_hat(:,:,6));
W_hat  = cat(3, S1.W_hat(:,:,1),  S4.W_hat(:,:,4),  S6.W_hat(:,:,6));

U2_hat = cat(3, S1.U2_hat(:,:,1), S4.U2_hat(:,:,4), S6.U2_hat(:,:,6));
V2_hat = cat(3, S1.V2_hat(:,:,1), S4.V2_hat(:,:,4), S6.V2_hat(:,:,6));
W2_hat = cat(3, S1.W2_hat(:,:,1), S4.W2_hat(:,:,4), S6.W2_hat(:,:,6));

% optional: sigma from each file
sigma_list = [S1.result_sigma(1,kz_idx), S4.result_sigma(1,kz_idx), S6.result_sigma(1,kz_idx)];

show_velocity = "true";
compute_singular_value = "false"
mesh_independence = "false"
c_number = 12; % 96 % MULTIPLE OF NUMBER OF CORES

% ---- choose the kz index you ran (adjust if needed) ----
                    % <--- set to the i you used in the run

load data_x_Re300.mat
load data_y_Re300.mat
load u_mean_zt_Re300.mat               % fine-grid mean U from Dedalus
load x_solid_Re300.mat
load y_solid_Re300.mat
filename_mask = ['mask_smooth_hillperiodic_Re' num2str(Re) '_' num2str(Ny) 'x' num2str(Nx) '.mat'];
load(filename_mask)
%load mask_smooth_hillperiodic_Re100_120x96.mat

solid_mask = (mask_smooth == 300);   % <<< IMPORTANT: > 0, not == 0
data_y = data_y + 1.5
y1_vals = 3*C * ( 1 + tanh( B*(abs(data_x - A) - B) ) );   % bottom wall y'(x,0)
y1_vals = y1_vals - 1.5;

[X_old,Y_old] = meshgrid(data_x, data_y);
[Xc, Yc]      = meshgrid(x, y);  % 'x','y' came from stability_results.mat

Umean_coarse  = interp2(X_old, Y_old, u_mean_zt, Xc, Yc, 'linear');  % note transpose
[dUdx_coarse, dUdy_coarse] = gradient(Umean_coarse, x, y);

% (optional) clear big vars:
clear u_mean_zt X_old Y_old Xc Yc


% make subfolder
snapdir = fullfile(folderpath, sprintf('snapshots_hillp_Re%d_kz%0.4g_multiC',Re, kz));
if ~exist(snapdir,'dir'), mkdir(snapdir); end


% grid for quiver (same grid you used for pcolor/contour)
[Xg, Yg] = meshgrid(x, y);

% thin the arrows a bit so the plot stays readable
qstep = max(1, round(min(numel(x), numel(y))/40));   % ~25–35 arrows each way
Xs = Xg(1:qstep:end, 1:qstep:end);
Ys = Yg(1:qstep:end, 1:qstep:end);

%guards for walls & colormap
haveWalls = exist('x_val','var') && exist('y1_vals','var') && exist('y2_vals','var');
hasBWR   = exist('bluewhitered','file');

% ---- layout + tick style (same idea as working case) ----
BIG_TICKS = 40;                         % axis tick fontsize
NXT       = 4;                          % number of major x ticks
NYT       = 4;                          % number of major y ticks

% FIG_POS   = [100 100 1350 1000];        % figure size in pixels
% AX_POS = [0.18 0.26 0.55 0.62]; % axes position [left bottom width height]
% CB_POS    = [0.78 0.16 0.05 0.78];      % colorbar position
% Lx        = 9.0;                        % for last x–tick

FIG_POS = [100 100 1350 1000];
AX_POS  = [0.12 0.28 0.78 0.58];
CB_POS = [0.18 0.20 0.64 0.04];
Lx        = 9.0;                        % for last x–tick


if show_velocity == "true"

    idx_list = 1:nc;

    Umax_common  = 0; U2max_common = 0;
    Vmax_common  = 0; V2max_common = 0;
    Wmax_common  = 0; W2max_common = 0;

    for ii = idx_list
        tmp = U_hat(:,:,ii);   tmp(solid_mask) = NaN;
        Umax_common = max(Umax_common, max(abs(tmp(:)), [], 'omitnan'));

        tmp = U2_hat(:,:,ii);  tmp(solid_mask) = NaN;
        U2max_common = max(U2max_common, max(abs(tmp(:)), [], 'omitnan'));

        tmp = V_hat(:,:,ii);   tmp(solid_mask) = NaN;
        Vmax_common = max(Vmax_common, max(abs(tmp(:)), [], 'omitnan'));

        tmp = V2_hat(:,:,ii);  tmp(solid_mask) = NaN;
        V2max_common = max(V2max_common, max(abs(tmp(:)), [], 'omitnan'));

        tmp = W_hat(:,:,ii);   tmp(solid_mask) = NaN;
        Wmax_common = max(Wmax_common, max(abs(tmp(:)), [], 'omitnan'));

        tmp = W2_hat(:,:,ii);  tmp(solid_mask) = NaN;
        W2max_common = max(W2max_common, max(abs(tmp(:)), [], 'omitnan'));
    end

    %% ========================================================
    %  LOOP OVER ALL c CASES AND EXPORT PLOTS
    % =========================================================
    for i = 1:nc

        c_label = c_index(i);

        % ---------- U response ----------
        fU1 = figure('Visible','off','Position',[100 100 1000 800]);
        U_plot = U_hat(:,:,i);
        U_plot(solid_mask) = NaN;

        contourf(x, y, U_plot, 40, 'LineWidth', 0.5);
        hold on
        shading interp
        colormap(bluewhitered);
        caxis([-Umax_common Umax_common]);

        %writematrix(U_plot, fullfile(snapdir, sprintf('U_response_c%02d_kz%g.csv', c_label, kz)));

        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);
        hold on

        hillp_tick_function;
        daspect([1 1 1])

        exportgraphics(fU1, fullfile(snapdir, sprintf('U_response_Re%d_c%02d_kz%g.png',Re, c_label, kz)), ...
            'Resolution', 600, 'BackgroundColor', 'white');
        close(fU1);

        % ---------- X forcing ----------
        fU3 = figure('Visible','off','Position',[100 100 1000 800]);
        U2_plot = U2_hat(:,:,i);
        U2_plot(solid_mask) = NaN;

        contourf(x, y, U2_plot, 40, 'LineWidth', 0.5);
        hold on
        shading interp
        colormap(bluewhitered);
        caxis([-U2max_common U2max_common]);

        %writematrix(U2_plot, fullfile(snapdir, sprintf('X_forcing_c%02d_kz%g.csv', c_label, kz)));

        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);
        hold on

        hillp_tick_function;
        daspect([1 1 1])

        exportgraphics(fU3, fullfile(snapdir, sprintf('X_forcing_Re%d_c%02d_kz%g.png',Re, c_label, kz)), ...
            'Resolution', 600, 'BackgroundColor', 'white');
        close(fU3);

        % ---------- V response ----------
        fV1 = figure('Visible','off','Position',[100 100 1000 800]);
        V_plot = V_hat(:,:,i);
        V_plot(solid_mask) = NaN;

        contourf(x, y, V_plot, 40, 'LineWidth', 0.5);
        hold on
        shading interp
        colormap(bluewhitered);
        caxis([-Vmax_common Vmax_common]);

        %writematrix(V_plot, fullfile(snapdir, sprintf('V_response_c%02d_kz%g.csv', c_label, kz)));

        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);
        hold on

        hillp_tick_function;
        daspect([1 1 1])

        exportgraphics(fV1, fullfile(snapdir, sprintf('V_response_Re%d_c%02d_kz%g.png',Re, c_label, kz)), ...
            'Resolution', 600, 'BackgroundColor', 'white');
        close(fV1);

        % ---------- Y forcing ----------
        fV3 = figure('Visible','off','Position',[100 100 1000 800]);
        V2_plot = V2_hat(:,:,i);
        V2_plot(solid_mask) = NaN;

        contourf(x, y, V2_plot, 40, 'LineWidth', 0.5);
        hold on
        shading interp
        colormap(bluewhitered);
        caxis([-V2max_common V2max_common]);

        %writematrix(V2_plot, fullfile(snapdir, sprintf('Y_forcing_c%02d_kz%g.csv', c_label, kz)));

        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);
        hold on

        hillp_tick_function;
        daspect([1 1 1])

        exportgraphics(fV3, fullfile(snapdir, sprintf('Y_forcing_Re%d_c%02d_kz%g.png',Re, c_label, kz)), ...
            'Resolution', 600, 'BackgroundColor', 'white');
        close(fV3);

        % ---------- W response ----------
        fW1 = figure('Visible','off','Position',[100 100 1000 800]);
        W_plot = W_hat(:,:,i);
        W_plot(solid_mask) = NaN;

        contourf(x, y, W_plot, 40, 'LineWidth', 0.5);
        hold on
        shading interp
        colormap(bluewhitered);
        caxis([-Wmax_common Wmax_common]);

        %writematrix(W_plot, fullfile(snapdir, sprintf('W_response_c%02d_kz%g.csv', c_label, kz)));

        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);
        hold on

        hillp_tick_function;
        daspect([1 1 1])

        exportgraphics(fW1, fullfile(snapdir, sprintf('W_response_Re%d_c%02d_kz%g.png',Re, c_label, kz)), ...
            'Resolution', 600, 'BackgroundColor', 'white');
        close(fW1);

        % ---------- Z forcing ----------
        fW3 = figure('Visible','off','Position',[100 100 1000 800]);
        W2_plot = W2_hat(:,:,i);
        W2_plot(solid_mask) = NaN;

        contourf(x, y, W2_plot, 40, 'LineWidth', 0.5);
        hold on
        shading interp
        colormap(bluewhitered);
        caxis([-W2max_common W2max_common]);

        %writematrix(W2_plot, fullfile(snapdir, sprintf('Z_forcing_c%02d_kz%g.csv', c_label, kz)));

        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);
        hold on

        hillp_tick_function;
        daspect([1 1 1])

        exportgraphics(fW3, fullfile(snapdir, sprintf('Z_forcing_Re%d_c%02d_kz%g.png',Re, c_label, kz)), ...
            'Resolution', 300, 'BackgroundColor', 'white');
        close(fW3);
    end
end

