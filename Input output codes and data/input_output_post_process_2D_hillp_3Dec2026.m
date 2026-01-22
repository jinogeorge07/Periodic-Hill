% Input Output Post Processing for 2 way coupled flow
clc; clear all; close all;
%% Add path for color profile---------------------------------------------------------
%folderpath = 'E:/dedalus/Wavy Wall/dedalus_local_Re75_3/snapshots_channel/mean_v_Re75.00_c1.00'
folderpath = 'E:/dedalus/Wavy Wall/flexible signed distance function/opt_erf_peclet/dedalus_local_Re190_nek5000_hill_ubulk_2D/snapshots_channel/mean_v_Re190.00_c1.00'
cd(folderpath)
addpath('/gpfs/homefs1/jig23007/ray_tracing_matlab/cbrewer2-master');

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
compute_singular_value = "false"
mesh_independence = "false"
c_number = 12; % 96 % MULTIPLE OF NUMBER OF CORES
% ---- choose the kz index you ran (adjust if needed) ----
% <--- set to the i you used in the run
kz = kz_list(kz_idx);

load data_x.mat
load data_y.mat
load u_mean_zt.mat               % fine-grid mean U from Dedalus
load x_solid.mat
load y_solid.mat
filename_mask = ['mask_smooth_hillperiodic_Re' num2str(Re) '_' num2str(Ny) 'x' num2str(Nx) '.mat'];
load(filename_mask)
%load mask_smooth_hillperiodic_Re100_120x96.mat

solid_mask = (mask_smooth == 100);   % <<< IMPORTANT: > 0, not == 0
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

if compute_singular_value == "true"

    Ix = speye(Nx); Iy = speye(Ny);
    I  = kron(Ix,Iy);  % sparse because Ix, Iy are sparse
    [~,w] = clencurt(Ny-1);
    w = w*((2+2*epsilon)/2); % normalize weighted matrix

    % w = spdiags(reshape(w,[],1),0,Ny,Ny);
    % w = kron(Ix, diag(w));     % weighted matrix
    w_sqrt = sqrt(w);
    w_sqrt = spdiags(w_sqrt',0,Ny,Ny);
    w_sqrt = kron(Ix, w_sqrt);

    w_sqrt_inv = 1./(sqrt(w)); % elementwise sqrt
    w_sqrt_inv = spdiags(w_sqrt_inv',0,Ny,Ny);
    w_sqrt_inv = kron(Ix, w_sqrt_inv);

    for c_index = 1:c_number
        vec = result_U_svd{c_index, kz_idx};  assert(numel(vec) == 3*N); % divide by weighting D matrix
        vec2= result_V_svd{c_index, kz_idx};  assert(numel(vec2)== 3*N);

        % RESPONSE
        U_hat(:,:,c_index)  = reshape(w_sqrt_inv*abs(vec(1:N)),            Ny, Nx);
        V_hat(:,:,c_index)  = reshape(w_sqrt_inv*abs(vec(N+1:2*N)),        Ny, Nx);
        W_hat(:,:,c_index)  = reshape(w_sqrt_inv*abs(vec(2*N+1:3*N)),      Ny, Nx);

        % FORCING
        U2_hat(:,:,c_index) = reshape(w_sqrt_inv*abs(vec2(1:N)),           Ny, Nx);
        V2_hat(:,:,c_index) = reshape(w_sqrt_inv*abs(vec2(N+1:2*N)),       Ny, Nx);
        W2_hat(:,:,c_index) = reshape(w_sqrt_inv*abs(vec2(2*N+1:3*N)),     Ny, Nx);
    end


    if kz_idx == 5 && mesh_independence == "true"

        if Ny == 120
            % Maximum Response mode value 120x96
            max_Uhat_120x96 = squeeze(max(max(U_hat, [], 1), [], 2));
            max_Vhat_120x96 = squeeze(max(max(V_hat, [], 1), [], 2));

            save('max_Uhat_120x96.mat','max_Uhat_120x96');
            save('max_Vhat_120x96.mat','max_Vhat_120x96');

        elseif Ny == 96
            % Maximum Response mode value
            max_Uhat_96x80 = squeeze(max(max(U_hat, [], 1), [], 2));
            max_Vhat_96x80 = squeeze(max(max(V_hat, [], 1), [], 2));
            save('max_Uhat_96x80.mat','max_Uhat_96x80');
            save('max_Vhat_96x80.mat','max_Vhat_96x80');

        elseif Ny == 60
            % Maximum Response mode value
            max_Uhat_60x48 = squeeze(max(max(U_hat, [], 1), [], 2));
            max_Vhat_60x48 = squeeze(max(max(V_hat, [], 1), [], 2));
            save('max_Uhat_60x48.mat','max_Uhat_60x48');
            save('max_Vhat_60x48.mat','max_Vhat_60x48');

        elseif Ny == 72
            % Maximum Response mode value
            max_Uhat_72x60 = squeeze(max(max(U_hat, [], 1), [], 2));
            max_Vhat_72x60 = squeeze(max(max(V_hat, [], 1), [], 2));
            save('max_Uhat_72x60.mat','max_Uhat_72x60');
            save('max_Vhat_72x60.mat','max_Vhat_72x60');

        elseif Ny == 132
            % Maximum Response mode value
            max_Uhat_132x108 = squeeze(max(max(U_hat, [], 1), [], 2));
            max_Vhat_132x108 = squeeze(max(max(V_hat, [], 1), [], 2));
            save('max_Uhat_132x108.mat','max_Uhat_132x108');
            save('max_Vhat_132x108.mat','max_Vhat_132x108');

        end

    end


end

% make subfolder
snapdir = fullfile(folderpath, sprintf('snapshots_hillp_kz%0.4g', kz));
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

FIG_POS   = [100 100 1350 1000];        % figure size in pixels
AX_POS = [0.18 0.26 0.55 0.62]; % axes position [left bottom width height]
CB_POS    = [0.78 0.16 0.05 0.78];      % colorbar position
Lx        = 9.0;                        % for last x–tick

if show_velocity == "true"

    for i = 1:min(12, size(U_hat,3))

        crit_contour = Umean_coarse - c_list(i);   % Ny x Nx
        % ---------- U: contour Response mode ----------
        fU1 = figure('Visible','off','Position',[100 100 1000 800]);

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


        % fill solid patch
        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        hold on

        hillp_tick_function;

        saveas(fU1, fullfile(snapdir, sprintf('U_response_c%02d_kz%g.png', i, kz)));
        % or exportgraphics if you prefer
        % exportgraphics(fU1, fullfile(snapdir, ...), 'Resolution', 300);

        close(fU1);

        % ---------- U: contour forcing mode ----------
        fU3 = figure('Visible','off','Position',[100 100 1000 800]);
        %contourf(x, y, U2_hat(:,:,i), 30, 'LineWidth', 1/2);
        % Copy U and blank out the solid region
        U2_plot = U2_hat(:,:,i);
        U2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        contourf(x, y, U2_plot, 40, 'LineWidth', 0.5);
        shading interp
        colormap(bluewhitered);                     % ensures the center color is white
        caxis([-max(abs(U2_plot(:))) max(abs(U2_plot(:)))]);  % ensures 0 is centered
        hold on

        % fill solid patch
        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        hold on

        hillp_tick_function;
        saveas(fU3, fullfile(snapdir, sprintf('X_forcing_c%02d_kz%g.png', i, kz)));

        %     exportgraphics(gcf, fullfile(snapdir, sprintf('X_forcing_c%02d_kz%g.png', i, kz)), 'Resolution', 300);
        %     close(fU3);

        % ---------- V: contour response mode----------
        fV1 = figure('Visible','off','Position',[100 100 1000 800]);
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

        hillp_tick_function;
        saveas(fV1, fullfile(snapdir, sprintf('V_response_c%02d_kz%g.png', i, kz)));

        % exportgraphics(gcf, fullfile(snapdir, sprintf('V_response_c%02d_kz%g.png', i, kz)), 'Resolution', 300);
        close(fV1);

        % ---------- V: contour forcing mode----------
        fV3 = figure('Visible','off','Position',[100 100 1000 800]);
        %contourf(x, y, V2_hat(:,:,i), 30, 'LineWidth', 1/2);

        %% Copy U and blank out the solid region
        V2_plot = V2_hat(:,:,i);
        V2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        contourf(x, y, V2_plot, 40, 'LineWidth', 0.5);
        shading interp
        colormap(bluewhitered);                     % ensures the center color is white
        caxis([-max(abs(V2_plot(:))) max(abs(V2_plot(:)))]);  % ensures 0 is centered
        hold on

        % fill solid patch
        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        hold on

        hillp_tick_function;
        saveas(fV3, fullfile(snapdir, sprintf('Y_forcing_c%02d_kz%g.png', i, kz)));

        %exportgraphics(gcf, fullfile(snapdir, sprintf('Y_forcing_c%02d_kz%g.png', i, kz)), 'Resolution', 300);
        close(fV3);


        % ---------- W: contour Response mode ----------
        fW1 = figure('Visible','off','Position',[100 100 1000 800]);
        %contourf(x, y, W_hat(:,:,i), 30, 'LineWidth', 1/2);

        %% Copy U and blank out the solid region
        W_plot = W_hat(:,:,i);
        W_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        contourf(x, y, W_plot, 40, 'LineWidth', 1/2);
        shading interp
        colormap(bluewhitered);                     % ensures the center color is white
        caxis([-max(abs(W_plot(:))) max(abs(W_plot(:)))]);  % ensures 0 is centered
        hold on

        % fill solid patch
        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        hold on

        hillp_tick_function;
        saveas(fW1, fullfile(snapdir, sprintf('W_response_c%02d_kz%g.png', i, kz)));
        %exportgraphics(gcf, fullfile(snapdir, sprintf('W_response_c%02d_kz%g.png', i, kz)), 'Resolution', 300);
        close(fW1);

        % ---------- Z: contour Forcing mode ----------
        fW3 = figure('Visible','off','Position',[100 100 1000 800]);
        %contourf(x, y, W2_hat(:,:,i), 30, 'LineWidth', 1/2);
        %% Copy U and blank out the solid region
        W2_plot = W2_hat(:,:,i);
        W2_plot(solid_mask) = NaN;      % hide solid, keep fluid only
        contourf(x, y, W2_plot, 40, 'LineWidth', 1/2);
        shading interp
        colormap(bluewhitered);                     % ensures the center color is white
        caxis([-max(abs(W2_plot(:))) max(abs(W2_plot(:)))]);  % ensures 0 is centered
        hold on

        % fill solid patch
        fill_patch_hillp;
        plot(data_x, y1_vals, 'k', 'LineWidth', 2);   % thick black curve
        hold on

        hillp_tick_function;
        saveas(fW3, fullfile(snapdir, sprintf('Z_forcing_c%02d_kz%g.png', i, kz)));
        %exportgraphics(gcf, fullfile(snapdir, sprintf('Z_forcing_c%02d_kz%g.png', i, kz)), 'Resolution', 300);
        close(fW3);


        %     if quiver == "ON"
        %       % ---------- (U,V) quiver ----------
        %     fQ = figure('Visible','off','Position',[100 100 1000 800]);
        %     % optional light background: magnitude of (u,v) for context
        %     contourf(x, y, hypot(U_hat(:,:,i), V_hat(:,:,i)), 20, 'LineStyle','none');
        %     colorbar; hold on
        %
        %     % quiver of (u,v) on downsampled grid
        %     quiver(Xs, Ys, ...
        %            U_hat(1:qstep:end, 1:qstep:end, i), ...
        %            V_hat(1:qstep:end, 1:qstep:end, i), ...
        %            'k');  % black arrows
        %     if hasBWR, colormap(bluewhitered); end
        %     xlabel('x','FontSize',16); ylabel('y','FontSize',16);
        %         title(sprintf('Quiver of (u,v), c\\_index=%d, kz=%g', i, kz));
        %
        %     % overlay your geometry markers if present
        %     if exist('x_solid','var') && exist('y_solid','var')
        %         plot(x_solid, y_solid, '.', 'Color', 'k');
        %     end
        %
        %     if haveWalls
        %         plot(x_val, y1_vals, 'k');
        %         plot(x_val, y2_vals, 'k');
        %     end
        %
        %     saveas(fQ, fullfile(snapdir, sprintf('UV_quiver_c%02d_kz%g.png', i, kz)));
        %     close(fQ);
        %     end

    end

end
sigma_list = result_sigma(:,kz_idx)
filename = ['sigma_list_Re' num2str(Re) '_kz' num2str(kz_idx) '.mat'];
save(filename, 'sigma_list');
