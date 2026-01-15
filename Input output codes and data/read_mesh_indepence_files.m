% Input Output Post Processing for 2 way coupled flow
clc; clear all; close all;
%% Add path for color profile---------------------------------------------------------
%folderpath = 'E:/dedalus/Wavy Wall/dedalus_local_Re75_3/snapshots_channel/mean_v_Re75.00_c1.00'  
folderpath = './'

%% Periodic hill and Domain Properties
%Ny = 132; Nx = 108; Nz = 2; % Nz to Nx for streamwise terms
%N = Ny*Nx;
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
rwb1 = bluewhitered; 
% figure;
% plot(x_val,y1_vals)
% hold on
% quiver = "OFF";
kz_idx = 5; 

U_48x60 = readmatrix('writefiles/U_hat_48x60.csv');
V_48x60 = readmatrix('writefiles/V_hat_48x60.csv');
W_48x60 = readmatrix('writefiles/W_hat_48x60.csv');

U_60x72 = readmatrix('writefiles/U_hat_60x72.csv');
V_60x72 = readmatrix('writefiles/V_hat_60x72.csv');
W_60x72 = readmatrix('writefiles/W_hat_60x72.csv');

U_80x96 = readmatrix('writefiles/U_hat_80x96.csv');
V_80x96 = readmatrix('writefiles/V_hat_80x96.csv');
W_80x96 = readmatrix('writefiles/W_hat_80x96.csv');

U_96x120 = readmatrix('writefiles/U_hat_96x120.csv');
V_96x120 = readmatrix('writefiles/V_hat_96x120.csv');
W_96x120 = readmatrix('writefiles/W_hat_96x120.csv');

U_108x132 = readmatrix('writefiles/U_hat_108x132.csv');
V_108x132 = readmatrix('writefiles/V_hat_108x132.csv');
W_108x132 = readmatrix('writefiles/W_hat_108x132.csv');

x_48x60 = readmatrix('writefiles/x_48x60.csv');
y_48x60 = readmatrix('writefiles/y_48x60.csv');

x_60x72 = readmatrix('writefiles/x_60x72.csv');
y_60x72 = readmatrix('writefiles/y_60x72.csv');

x_80x96 = readmatrix('writefiles/x_80x96.csv');
y_80x96 = readmatrix('writefiles/y_80x96.csv');

x_96x120 = readmatrix('writefiles/x_96x120.csv');
y_96x120 = readmatrix('writefiles/y_96x120.csv');

x_108x132 = readmatrix('writefiles/x_108x132.csv');
y_108x132 = readmatrix('writefiles/y_108x132.csv');

% 48 × 60
[X_48x60, Y_48x60] = meshgrid(x_48x60, y_48x60);

% 60 × 72
[X_60x72, Y_60x72] = meshgrid(x_60x72, y_60x72);

% 80 × 96
[X_80x96, Y_80x96] = meshgrid(x_80x96, y_80x96);

% 96 × 120
[X_96x120, Y_96x120] = meshgrid(x_96x120, y_96x120);

% 108 × 132
[X_108x132, Y_108x132] = meshgrid(x_108x132, y_108x132);


% 108x132 -> 48x60 
U_108to48x60 = interp2(X_108x132, Y_108x132, U_108x132, X_48x60, Y_48x60, 'linear', NaN);
V_108to48x60 = interp2(X_108x132, Y_108x132, V_108x132, X_48x60, Y_48x60, 'linear', NaN);
W_108to48x60 = interp2(X_108x132, Y_108x132, W_108x132, X_48x60, Y_48x60, 'linear', NaN);

% 108x132 -> 60x72
U_108to60x72 = interp2(X_108x132, Y_108x132, U_108x132, X_60x72, Y_60x72, 'linear', NaN);
V_108to60x72 = interp2(X_108x132, Y_108x132, V_108x132, X_60x72, Y_60x72, 'linear', NaN);
W_108to60x72 = interp2(X_108x132, Y_108x132, W_108x132, X_60x72, Y_60x72, 'linear', NaN);

% 108x132 -> 80x96
U_108to80x96 = interp2(X_108x132, Y_108x132, U_108x132, X_80x96, Y_80x96, 'linear', NaN);
V_108to80x96 = interp2(X_108x132, Y_108x132, V_108x132, X_80x96, Y_80x96, 'linear', NaN);
W_108to80x96 = interp2(X_108x132, Y_108x132, W_108x132, X_80x96, Y_80x96, 'linear', NaN);

% 108x132 -> 96x120
U_108to96x120 = interp2(X_108x132, Y_108x132, U_108x132, X_96x120, Y_96x120, 'linear', NaN);
V_108to96x120 = interp2(X_108x132, Y_108x132, V_108x132, X_96x120, Y_96x120, 'linear', NaN);
W_108to96x120 = interp2(X_108x132, Y_108x132, W_108x132, X_96x120, Y_96x120, 'linear', NaN);

den_U_108 = sqrt(mean(U_108x132(:).^2, 'omitnan'));
den_V_108 = sqrt(mean(V_108x132(:).^2, 'omitnan'));
den_W_108 = sqrt(mean(W_108x132(:).^2, 'omitnan'));


%% 108 -> 48
mask_48 = isfinite(U_108to48x60) & isfinite(V_108to48x60) & isfinite(W_108to48x60) & ...
          isfinite(U_48x60)      & isfinite(V_48x60)      & isfinite(W_48x60);

err_U_48 = sqrt(mean((U_48x60(mask_48) - U_108to48x60(mask_48)).^2)) / den_U_108;
err_V_48 = sqrt(mean((V_48x60(mask_48) - V_108to48x60(mask_48)).^2)) / den_V_108;
err_W_48 = sqrt(mean((W_48x60(mask_48) - W_108to48x60(mask_48)).^2)) / den_W_108;


%% 108 -> 60
mask_60 = isfinite(U_108to60x72) & isfinite(V_108to60x72) & isfinite(W_108to60x72) & ...
          isfinite(U_60x72)      & isfinite(V_60x72)      & isfinite(W_60x72);

err_U_60 = sqrt(mean((U_60x72(mask_60) - U_108to60x72(mask_60)).^2)) / den_U_108;
err_V_60 = sqrt(mean((V_60x72(mask_60) - V_108to60x72(mask_60)).^2)) / den_V_108;
err_W_60 = sqrt(mean((W_60x72(mask_60) - W_108to60x72(mask_60)).^2)) / den_W_108;

%% 108 -> 80
mask_80 = isfinite(U_108to80x96) & isfinite(V_108to80x96) & isfinite(W_108to80x96) & ...
          isfinite(U_80x96)      & isfinite(V_80x96)      & isfinite(W_80x96);

err_U_80 = sqrt(mean((U_80x96(mask_80) - U_108to80x96(mask_80)).^2)) / den_U_108;
err_V_80 = sqrt(mean((V_80x96(mask_80) - V_108to80x96(mask_80)).^2)) / den_V_108;
err_W_80 = sqrt(mean((W_80x96(mask_80) - W_108to80x96(mask_80)).^2)) / den_W_108;


%% 108 -> 96
mask_96 = isfinite(U_108to96x120) & isfinite(V_108to96x120) & isfinite(W_108to96x120) & ...
          isfinite(U_96x120)      & isfinite(V_96x120)      & isfinite(W_96x120);

err_U_96 = sqrt(mean((U_96x120(mask_96) - U_108to96x120(mask_96)).^2)) / den_U_108;
err_V_96 = sqrt(mean((V_96x120(mask_96) - V_108to96x120(mask_96)).^2)) / den_V_108;
err_W_96 = sqrt(mean((W_96x120(mask_96) - W_108to96x120(mask_96)).^2)) / den_W_108;


fprintf('Errors vs 108 (denom = native 108):\n');
fprintf('48 : U %.3e  V %.3e  W %.3e\n', err_U_48, err_V_48, err_W_48);
fprintf('60 : U %.3e  V %.3e  W %.3e\n', err_U_60, err_V_60, err_W_60);
fprintf('80 : U %.3e  V %.3e  W %.3e\n', err_U_80, err_V_80, err_W_80);
fprintf('96 : U %.3e  V %.3e  W %.3e\n', err_U_96, err_V_96, err_W_96);

% Reverse order: 108 -> 96 -> 80 -> 60
grids = [108 96 80 60];

errU = [err_U_96, err_U_80, err_U_60, err_U_48];
errV = [err_V_96, err_V_80, err_V_60, err_V_48];
errW = [err_W_96, err_W_80, err_W_60, err_W_48];


outdir = 'writefiles';   % or full path like 'C:/Users/you/Documents/results/'
if ~exist(outdir, 'dir')
    mkdir(outdir);
end


figure(49)
plot(grids, errU, '-o', 'LineWidth', 2); hold on

xlabel('Fine grid resolution');
ylabel('Relative RMS error vs 48x60');
title('Grid convergence (fine \rightarrow coarse)');
legend('U','Location','northeast');
grid on
saveas(gcf, fullfile(outdir,'U_grid_convergence_linear_fine_to_coarse.png'));


figure(50)
plot(grids, errU, '-o', 'LineWidth', 2); hold on
plot(grids, errV, '-s', 'LineWidth', 2);hold on
plot(grids, errW, '-^', 'LineWidth', 2); hold on

xlabel('Fine grid resolution');
ylabel('Relative RMS error vs 48x60');
title('Grid convergence (fine \rightarrow coarse)');
legend('U','V','W','Location','northeast');
grid on
saveas(gcf, fullfile(outdir,'UVW_grid_convergence_linear_fine_to_coarse.png'));

