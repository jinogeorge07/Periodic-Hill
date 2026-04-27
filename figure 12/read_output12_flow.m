% Input Output Post Processing for 2 way coupled flow
clc; clear all; close all;
%% Add path for color profile---------------------------------------------------------
%=========================================================================================================================
folderpath = './'

%% Periodic hill and Domain Properties
Ny = 132; Nx = 108; Nz = 2; % Nz to Nx for streamwise terms
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
Re = 190
x_val = linspace(0, Lx, Nx);  % replaces fourdif x values for plotting only
A  = 4.5;
B  = 3.5;
C  = 1/6;

% figure;
% plot(x_val,y1_vals)
% hold on
% quiver = "OFF";
kz_idx = 20;
load('stability_results_Re190_132x108_kz20.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega

%load('stability_results_kz32_60x48.mat');
%load('stability_results.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega

show_velocity = "true";
c_number = 12; % 96 % MULTIPLE OF NUMBER OF CORES
% ---- choose the kz index you ran (adjust if needed) ----
% <--- set to the i you used in the run
kz = kz_list(kz_idx);
kz = 0.55;

load data_x.mat
load data_y.mat
load u_mean_zt.mat               % fine-grid mean U from Dedalus
load x_solid.mat
load y_solid.mat
filename_mask = ['mask_smooth_hillperiodic_Re' num2str(Re) '_' num2str(Ny) 'x' num2str(Nx) '.mat'];
load(filename_mask)
rwb1 = bluewhitered
%load mask_smooth_hillperiodic_Re100_120x96.mat

solid_mask = (mask_smooth == Re);   % <<< IMPORTANT: > 0, not == 0
data_y = data_y + 1.5
y1_vals = 3*C * ( 1 + tanh( B*(abs(data_x - A) - B) ) );   % bottom wall y'(x,0)
y1_vals = y1_vals - 1.5;

[X_old,Y_old] = meshgrid(data_x, data_y);
[Xc, Yc]      = meshgrid(x, y);  % 'x','y' came from stability_results.mat

Umean_coarse  = interp2(X_old, Y_old, u_mean_zt, Xc, Yc, 'linear');  % note transpose
[dUdx_coarse, dUdy_coarse] = gradient(Umean_coarse, x, y);

% Shear layer
% define a "large shear" mask from |dUdy|; tweak 0.4 if needed
th_shear   = 0.4 * max(abs(dUdy_coarse(:)));      % 40% of max shear
shear_mask = abs(dUdy_coarse) > th_shear;

% (optional) clear big vars:
clear u_mean_zt X_old Y_old Xc Yc

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

FIG_POS   = [100 100 1350 1000];        % figure size in pixels
AX_POS  = [0.12 0.28 0.78 0.58];
CB_POS = [0.18 0.20 0.64 0.04];
Lx        = 9.0;                        % for last x–tick

if show_velocity == "true"

    idx_list = [1,4,6];

    Umax_common = 0;
    for ii = idx_list
        tmp = U2_hat(:,:,ii);
        tmp(solid_mask) = NaN;
        Umax_common = max(Umax_common, max(abs(tmp(:)), [], 'omitnan'));
    end
    Vmax_common = 0;
    for ii = idx_list
        tmp = V2_hat(:,:,ii);
        tmp(solid_mask) = NaN;
        Vmax_common = max(Vmax_common, max(abs(tmp(:)), [], 'omitnan'));
    end

    Wmax_common = 0;
    for ii = idx_list
        tmp = W2_hat(:,:,ii);
        tmp(solid_mask) = NaN;
        Wmax_common = max(Wmax_common, max(abs(tmp(:)), [], 'omitnan'));
    end

    for i = [1,4,6]

        %     % ---------- U: contour Response mode ----------
        %     fU1 = figure('Visible','off','Position',[100 100 1000 800]);
        %
        %     % Copy U and blank out the solid region
        %     U_plot = U_hat(:,:,i);
        %     U_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        %     contourf(x, y, U_plot, 40, 'LineWidth', 0.5);
        %     hold on
        %     shading interp
        %     colormap(bluewhitered);                     % ensures the center color is white
        %     %contour(x, y, crit_contour, [0 0], 'k', 'LineWidth', 2);   % critical layer
        %     caxis([-max(abs(U_plot(:))) max(abs(U_plot(:)))]);  % ensures 0 is centered
        %     hold on
        %
        %
        %     % fill solid patch
        %     fill_patch_hillp;
        %     plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        %     hold on
        %
        %     hillp_tick_function;
        %
        %     saveas(fU1, fullfile(snapdir, sprintf('U_response_c%02d_kz%g.png', i, kz)));
        %     % or exportgraphics if you prefer
        %     % exportgraphics(fU1, fullfile(snapdir, ...), 'Resolution', 300);
        %
        %     close(fU1);

        % ---------- X: contour forcing mode ----------
        fU3 = figure('Visible','on','Position',[100 100 1000 800]);
        %contourf(x, y, U2_hat(:,:,i), 30, 'LineWidth', 1/2);
        % Copy U and blank out the solid region
        U2_plot = U2_hat(:,:,i);
        U2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        contourf(x, y, U2_plot, 40, 'LineWidth', 0.5);
        shading interp
        colormap(bluewhitered);                     % ensures the center color is white
        %caxis([-max(abs(U2_plot(:))) max(abs(U2_plot(:)))]);  % ensures 0 is centered
        caxis([-Umax_common Umax_common]);
        hold on

        writematrix(U2_plot, fullfile(sprintf('X_forcing_c%02d_kz%g.csv', i, kz)));
        % fill solid patch
        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        hold on

        hillp_tick_function;
        daspect([1 1 1]);
        %saveas(fU3, fullfile(sprintf('X_forcing_c%02d_kz%g_Re%d.png', i, kz,Re)));
        exportgraphics(fU3, sprintf('X_forcing_c%02d_kz%g_Re%d.png', i, kz,Re), ...
        'Resolution', 600, 'BackgroundColor', 'white');       
        close(fU3);

%             % ---------- V: contour response mode----------
%             fV1 = figure('Visible','off','Position',[100 100 1000 800]);
%             %contourf(x, y, V_hat(:,:,i), 30, 'LineWidth', 1/2);
%         
%             %% Copy U and blank out the solid region
%             % Copy U and blank out the solid region
%             V_plot = V_hat(:,:,i);
%             V_plot(solid_mask) = NaN;      % hide solid, keep fluid only
%             contourf(x, y, V_plot, 40, 'LineWidth', 0.5);
%             shading interp
%             colormap(bluewhitered);                     % ensures the center color is white
%             %hold on
%             %contour(x, y, crit_contour, [0 0], 'k', 'LineWidth', 2);   % critical layer
%             caxis([-max(abs(V_plot(:))) max(abs(V_plot(:)))]);  % ensures 0 is centered
%             hold on
%         
%             % fill solid patch
%             fill_patch_hillp;
%             plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
%             hold on
%         
%             hillp_tick_function;
%             saveas(fV1, fullfile(snapdir, sprintf('V_response_c%02d_kz%g.png', i, kz)));
%         
%            % exportgraphics(gcf, fullfile(snapdir, sprintf('V_response_c%02d_kz%g.png', i, kz)), 'Resolution', 300);
        %     close(fV1);

        % ---------- Y: contour forcing mode----------
        fV3 = figure('Visible','on','Position',[100 100 1000 800]);
        %contourf(x, y, V2_hat(:,:,i), 30, 'LineWidth', 1/2);

        %% Copy U and blank out the solid region
        V2_plot = V2_hat(:,:,i);
        V2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        contourf(x, y, V2_plot, 40, 'LineWidth', 0.5);
        shading interp
        colormap(bluewhitered);                     % ensures the center color is white
        %caxis([-max(abs(V2_plot(:))) max(abs(V2_plot(:)))]);  % ensures 0 is centered
        caxis([-Vmax_common Vmax_common]);
        hold on

        writematrix(V2_plot, fullfile(sprintf('Y_forcing_c%02d_kz%g.csv', i, kz)));
        % fill solid patch
        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        hold on

        hillp_tick_function;
        daspect([1 1 1]);
        %saveas(fV3, fullfile(sprintf('Y_forcing_c%02d_kz%g_Re%d.png', i, kz,Re)));

        exportgraphics(fV3, sprintf('Y_forcing_c%02d_kz%g_Re%d.png', i, kz,Re), ...
        'Resolution', 600, 'BackgroundColor', 'white');              
        close(fV3);


%             % ---------- W: contour Response mode ----------
%             fW1 = figure('Visible','off','Position',[100 100 1000 800]);
%             %contourf(x, y, W_hat(:,:,i), 30, 'LineWidth', 1/2);
%         
%             %% Copy U and blank out the solid region
%             W_plot = W_hat(:,:,i);
%             W_plot(solid_mask) = NaN;      % hide solid, keep fluid only
%             contourf(x, y, W_plot, 40, 'LineWidth', 1/2);
%             shading interp
%             colormap(bluewhitered);                     % ensures the center color is white
%             caxis([-max(abs(W_plot(:))) max(abs(W_plot(:)))]);  % ensures 0 is centered
%             hold on
%         
%             % fill solid patch
%             fill_patch_hillp;
%             plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
%             hold on
%         
%             hillp_tick_function;
%             saveas(fW1, fullfile(snapdir, sprintf('W_response_c%02d_kz%g.png', i, kz)));
%             %exportgraphics(gcf, fullfile(snapdir, sprintf('W_response_c%02d_kz%g.png', i, kz)), 'Resolution', 300);
%             close(fW1);

%         % ---------- Z: contour Forcing mode ----------
%         fW3 = figure('Visible','on','Position',[100 100 1000 800]);
%         %contourf(x, y, W2_hat(:,:,i), 30, 'LineWidth', 1/2);
%         %% Copy U and blank out the solid region
%         W2_plot = W2_hat(:,:,i);
%         W2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
%         contourf(x, y, W2_plot, 40, 'LineWidth', 1/2);
%         shading interp
%         colormap(bluewhitered);                     % ensures the center color is white
%         %caxis([-max(abs(W2_plot(:))) max(abs(W2_plot(:)))]);  % ensures 0 is centered
%         caxis([-Wmax_common Wmax_common]);
%         hold on
% 
%         % fill solid patch
%         fill_patch_hillp;
%         plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
%         hold on
% 
%         hillp_tick_function;
%         daspect([1 1 1]);
%         %saveas(fW3, fullfile(sprintf('Z_forcing_c%02d_kz%g_Re%d.png', i, kz,Re)));
%         exportgraphics(fW3, sprintf('Z_forcing_c%02d_kz%g_Re%d.png', i, kz,Re), ...
%         'Resolution', 600, 'BackgroundColor', 'white');               
%         close(fW3);

    end

end

sigma_list = result_sigma(:,kz_idx)
filename = ['sigma_list_Re' num2str(Re) '_kz' num2str(kz_idx) '.mat'];
save(filename, 'sigma_list');

%% ------------------------------------------------------------------------
%% Re = 300 ---------------------------------------------------------------

Ny = 168; Nx = 120; Nz = 2; % Nz to Nx for streamwise terms
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
%kz = kz_list(kz_idx);
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

% % make subfolder
% snapdir = fullfile(folderpath, sprintf('snapshots_hillp_Re%d_kz%0.4g_multiC',Re, kz));
% if ~exist(snapdir,'dir'), mkdir(snapdir); end

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

        exportgraphics(fU3, fullfile(sprintf('X_forcing_c%02d_kz%g_Re%d.png', c_label, kz, Re)), ...
            'Resolution', 600, 'BackgroundColor', 'white');
        close(fU3);

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

        exportgraphics(fV3, fullfile(sprintf('Y_forcing_c%02d_kz%g_Re%d.png',c_label, kz, Re)), ...
            'Resolution', 600, 'BackgroundColor', 'white');
        close(fV3);

    end
end


