% Input Output Post Processing for 2 way coupled flow
clc; clear all; close all;
%% Add path for color profile---------------------------------------------------------
folderpath = './'

%% Periodic hill and Domain Properties
Ny = 120; Nx = 96; Nz = 2; % Nz to Nx for streamwise terms
N = Ny*Nx;
epsilon = 0.50;
%epsilon = 0.25;
h = 1;
y0 = h;
A1 = epsilon;
A2 = epsilon;
%Lx = 0.6 * pi;
Lx = 9.0;
Ly = 1.0 + 1.0*epsilon;
Re = 100
x_val = linspace(0, Lx, Nx);  % replaces fourdif x values for plotting only
A  = 4.5;
B  = 3.5;
C  = 1/6;

% figure;
% plot(x_val,y1_vals)
% hold on
% quiver = "OFF";
kz_idx = [5,20,32];
% load('stability_results_Re100_120x96_kz5.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega
% load('stability_results_Re100_120x96_kz20.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega
% load('stability_results_Re100_120x96_kz32.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega

%load('stability_results_kz32_60x48.mat');
%load('stability_results.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega

c_number = 12; % 96 % MULTIPLE OF NUMBER OF CORES
% ---- choose the kz index you ran (adjust if needed) ----
% <--- set to the i you used in the run

load data_x.mat
load data_y.mat
load u_mean_zt.mat               % fine-grid mean U from Dedalus
load x_solid.mat
load y_solid.mat
filename_mask = ['mask_smooth_hillperiodic_Re' num2str(Re) '_' num2str(Ny) 'x' num2str(Nx) '.mat'];
load(filename_mask)
rwb1 = bluewhitered
%load mask_smooth_hillperiodic_Re100_120x96.mat

solid_mask = (mask_smooth == 100);   % <<< IMPORTANT: > 0, not == 0
data_y = data_y + 1.5
y1_vals = 3*C * ( 1 + tanh( B*(abs(data_x - A) - B) ) );   % bottom wall y'(x,0)
y1_vals = y1_vals - 1.5;

%guards for walls & colormap
haveWalls = exist('x_val','var') && exist('y1_vals','var') && exist('y2_vals','var');
hasBWR   = exist('bluewhitered','file');

% ---- layout + tick style (same idea as working case) ----
BIG_TICKS = 40;                         % axis tick fontsize
NXT       = 4;                          % number of major x ticks
NYT       = 4;                          % number of major y ticks

FIG_POS   = [100 100 1350 1000];        % figure size in pixels
AX_POS = [0.18 0.26 0.55 0.62]; % axes position [left bottom width height]
CB_POS    = [0.78 0.16 0.05 0.78];      % colorbar position
Lx        = 9.0;                        % for last x–tick

for kz_idx = [5,20,32]

    fname = ['stability_results_Re100_120x96_kz' num2str(kz_idx) '.mat'];
    load(fname);

    kz = kz_list(kz_idx);

    [X_old,Y_old] = meshgrid(data_x, data_y);
    [Xc, Yc]      = meshgrid(x, y);  % 'x','y' came from stability_results.mat

    Umean_coarse  = interp2(X_old, Y_old, u_mean_zt, Xc, Yc, 'linear');  % note transpose
    [dUdx_coarse, dUdy_coarse] = gradient(Umean_coarse, x, y);

    % (optional) clear big vars:
    %clear u_mean_zt X_old Y_old Xc Yc

    % grid for quiver (same grid you used for pcolor/contour)
    [Xg, Yg] = meshgrid(x, y);

    % thin the arrows a bit so the plot stays readable
    qstep = max(1, round(min(numel(x), numel(y))/40));   % ~25–35 arrows each way
    Xs = Xg(1:qstep:end, 1:qstep:end);
    Ys = Yg(1:qstep:end, 1:qstep:end);

    for i = [1]

        crit_contour = Umean_coarse - c_list(i);   % Ny x Nx

        % ---------- U: contour Response mode ----------
        fU1 = figure('Visible','on','Position',[100 100 1000 800]);
        % Copy U and blank out the solid region
        U_plot = U_hat(:,:,i);
        U_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        contourf(x, y, U_plot, 40, 'LineWidth', 0.5);
        hold on
        shading interp
        colormap(bluewhitered);                     % ensures the center color is white
        %contour(x, y, crit_contour, [0 0], 'k', 'LineWidth', 2);   % critical layer
        caxis([-max(abs(U_plot(:))) max(abs(U_plot(:)))]);  % ensures 0 is centered
        hold on

        writematrix(U_plot, fullfile(sprintf('U_response_c%02d_kz%g.csv', i, kz)));
        % fill solid patch
        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        hold on

        hillp_tick_function;

        saveas(fU1, fullfile(sprintf('U_response_c%02d_kz%g.png', i, kz)));
        close(fU1);

        % %     % ---------- U: contour forcing mode ----------
        %     fU3 = figure('Visible','off','Position',[100 100 1000 800]);
        %     %contourf(x, y, U2_hat(:,:,i), 30, 'LineWidth', 1/2);
        %     % Copy U and blank out the solid region
        %     U2_plot = U2_hat(:,:,i);
        %     U2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        %     contourf(x, y, U2_plot, 40, 'LineWidth', 0.5);
        %     shading interp
        %     colormap(bluewhitered);                     % ensures the center color is white
        %     caxis([-max(abs(U2_plot(:))) max(abs(U2_plot(:)))]);  % ensures 0 is centered
        %     hold on
        %
        %     % fill solid patch
        %     fill_patch_hillp;
        %     plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        %     hold on
        %
        %     hillp_tick_function;
        %     saveas(fU3, fullfile(sprintf('X_forcing_c%02d_kz%g.png', i, kz)));
        %     close(fU3);


        % ---------- V: contour response mode----------
        fV1 = figure('Visible','on','Position',[100 100 1000 800]);
        %contourf(x, y, V_hat(:,:,i), 30, 'LineWidth', 1/2);

        %% Copy U and blank out the solid region
        % Copy U and blank out the solid region
        V_plot = V_hat(:,:,i);
        V_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        contourf(x, y, V_plot, 40, 'LineWidth', 0.5);
        shading interp
        colormap(bluewhitered);                     % ensures the center color is white
        %hold on
        %contour(x, y, crit_contour, [0 0], 'k', 'LineWidth', 2);   % critical layer
        caxis([-max(abs(V_plot(:))) max(abs(V_plot(:)))]);  % ensures 0 is centered
        hold on

        % fill solid patch
        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        hold on

        writematrix(V_plot, fullfile(sprintf('V_response_c%02d_kz%g.csv', i, kz)));
        hillp_tick_function;
        saveas(fV1, fullfile(sprintf('V_response_c%02d_kz%g.png', i, kz)));
        close(fV1);

        % %     % ---------- V: contour forcing mode----------
        %     fV3 = figure('Visible','off','Position',[100 100 1000 800]);
        %     %contourf(x, y, V2_hat(:,:,i), 30, 'LineWidth', 1/2);
        %
        %     %% Copy U and blank out the solid region
        %     V2_plot = V2_hat(:,:,i);
        %     V2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        %     contourf(x, y, V2_plot, 40, 'LineWidth', 0.5);
        %     shading interp
        %     colormap(bluewhitered);                     % ensures the center color is white
        %     caxis([-max(abs(V2_plot(:))) max(abs(V2_plot(:)))]);  % ensures 0 is centered
        %     hold on
        %
        %     % fill solid patch
        %     fill_patch_hillp;
        %     plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        %     hold on
        %
        %     hillp_tick_function;
        %     saveas(fV3, fullfile(sprintf('Y_forcing_c%02d_kz%g.png', i, kz)));
        %     close(fV3);


        %     % ---------- W: contour Response mode ----------
        %     fW1 = figure('Visible','off','Position',[100 100 1000 800]);
        %     %contourf(x, y, W_hat(:,:,i), 30, 'LineWidth', 1/2);
        %
        %     %% Copy U and blank out the solid region
        %     W_plot = W_hat(:,:,i);
        %     W_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        %     contourf(x, y, W_plot, 40, 'LineWidth', 1/2);
        %     shading interp
        %     colormap(bluewhitered);                     % ensures the center color is white
        %     caxis([-max(abs(W_plot(:))) max(abs(W_plot(:)))]);  % ensures 0 is centered
        %     hold on
        %
        %     % fill solid patch
        %     fill_patch_hillp;
        %     plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        %     hold on
        %
        %     hillp_tick_function;
        %     saveas(fW1, fullfile(sprintf('W_response_c%02d_kz%g.png', i, kz)));
        %     close(fW1);

        % %     % ---------- Z: contour Forcing mode ----------
        %     fW3 = figure('Visible','off','Position',[100 100 1000 800]);
        %     %contourf(x, y, W2_hat(:,:,i), 30, 'LineWidth', 1/2);
        %     %% Copy U and blank out the solid region
        %     W2_plot = W2_hat(:,:,i);
        %     W2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        %     contourf(x, y, W2_plot, 40, 'LineWidth', 1/2);
        %     shading interp
        %     colormap(bluewhitered);                     % ensures the center color is white
        %     caxis([-max(abs(W2_plot(:))) max(abs(W2_plot(:)))]);  % ensures 0 is centered
        %     hold on
        %
        %     % fill solid patch
        %     fill_patch_hillp;
        %     plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        %     hold on
        %
        %     hillp_tick_function;
        %     saveas(fW3, fullfile(sprintf('Z_forcing_c%02d_kz%g.png', i, kz)));
        %     close(fW3);

    end

end
% sigma_list = result_sigma(:,kz_idx)
% filename = ['sigma_list_Re' num2str(Re) '_kz' num2str(kz_idx) '.mat'];
% save(filename, 'sigma_list');
