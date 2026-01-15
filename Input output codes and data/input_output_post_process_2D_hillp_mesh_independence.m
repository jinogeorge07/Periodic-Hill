% Input Output Post Processing for 2 way coupled flow
clc; clear all; close all;
%% Add path for color profile---------------------------------------------------------
%folderpath = 'E:/dedalus/Wavy Wall/dedalus_local_Re75_3/snapshots_channel/mean_v_Re75.00_c1.00'  
folderpath = 'E:/dedalus/Wavy Wall/flexible signed distance function/opt_erf_peclet/dedalus_local_Re100_nek5000_hill_2a_ubulk_2D/snapshots_channel/mean_v_Re100.00_c1.00'
cd(folderpath)
addpath('/gpfs/homefs1/jig23007/ray_tracing_matlab/cbrewer2-master');

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

% figure;
% plot(x_val,y1_vals)
% hold on
% quiver = "OFF";
kz_idx = 5; 

outdir = 'writefiles';   % or full path like 'C:/Users/you/Documents/results/'
if ~exist(outdir, 'dir')
    mkdir(outdir);
end


for i=1:5

if i==5
    load('stability_results_Re100_132x108_kz5.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega
    Nx = 108; Ny = 132; N = Ny*Nx;
elseif i == 4
    load('stability_results_Re100_120x96_kz5.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega
    Nx = 96; Ny = 120; N = Ny*Nx;
elseif i == 3
    load('stability_results_Re100_96x80_kz5.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega
    Nx = 80; Ny = 96; N = Ny*Nx;
elseif i == 2
    load('stability_results_Re100_72x60_kz5.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega
    Nx = 60; Ny = 72; N = Ny*Nx;
elseif i == 1
    load('stability_results_Re100_60x48_kz5.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega
    Nx = 48; Ny = 60; N = Ny*Nx;

end

x_val = linspace(0, Lx, Nx);  % replaces fourdif x values for plotting only
A  = 4.5;
B  = 3.5;
C  = 1/6;

%load('stability_results_kz32_60x48.mat'); 
%load('stability_results.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega

show_velocity = "true";
compute_singular_value = "true"
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


% (optional) clear big vars:
clear u_mean_zt X_old Y_old Xc Yc


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

    for c_index = 6:6 % change from 1:c_number to omega = 0.55
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
    U_hat = squeeze(U_hat(:,:,6));V_hat = squeeze(V_hat(:,:,6));W_hat = squeeze(W_hat(:,:,6));
    U2_hat = squeeze(U2_hat(:,:,6));V2_hat = squeeze(V2_hat(:,:,6));W2_hat = squeeze(W2_hat(:,:,6));
   
    writematrix(U_hat,  fullfile(outdir, sprintf('U_hat_%dx%d.csv',  Nx, Ny)));
    writematrix(V_hat,  fullfile(outdir, sprintf('V_hat_%dx%d.csv',  Nx, Ny)));
    writematrix(W_hat,  fullfile(outdir, sprintf('W_hat_%dx%d.csv',  Nx, Ny)));

    writematrix(U2_hat, fullfile(outdir, sprintf('U2_hat_%dx%d.csv', Nx, Ny)));
    writematrix(V2_hat, fullfile(outdir, sprintf('V2_hat_%dx%d.csv', Nx, Ny)));
    writematrix(W2_hat, fullfile(outdir, sprintf('W2_hat_%dx%d.csv', Nx, Ny)));

    writematrix(x,  fullfile(outdir, sprintf('x_%dx%d.csv',  Nx, Ny)));
    writematrix(y,  fullfile(outdir, sprintf('y_%dx%d.csv',  Nx, Ny)));

end

sigma_list = result_sigma(:,kz_idx)
filename = ['sigma_list_Re' num2str(Re) '_kz' num2str(kz_idx) '.mat'];
save(filename, 'sigma_list');
