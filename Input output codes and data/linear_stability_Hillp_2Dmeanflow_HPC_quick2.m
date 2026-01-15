%% Code by Jino George . Group of Dr. Chang Liu University of Connecticut
%% same as linear_stability_wavy_wall_ChannelFlow_2Dmeanflow_HPC_June29_25,
%% with targeted performance edits (1–5) only.

clc; close all; clear all;

% --- Choose data folder by OS (minimal edits) ---
if isunix  % Linux/HPC
    folderpath = pwd;  % your job's scratch working directory
    addpath('/gpfs/homefs1/jig23007/ray_tracing_matlab/cbrewer2-master');
elseif ispc  % Windows/local
    %folderpath = 'E:\dedalus\Wavy Wall\flexible signed distance function\opt_erf_peclet\dedalus_local_wavy_wall_Re290_a1\snapshots_channel\mean_v_Re290.00_c1.00';
    folderpath = 'E:/dedalus/Wavy Wall/flexible signed distance function/opt_erf_peclet/dedalus_local_Re100_nek5000_hill_2a_ubulk_2D/snapshots_channel/mean_v_Re100.00_c1.00'
else
    error('Unsupported OS.');
end
if ~isfolder(folderpath), error('Data folder not found: %s', folderpath); end
cd(folderpath);

%% -------------CHOOSE SOLVER ---------------------------------------------
solver_type = 'svds';    % or 'svd'
solver_type_const = parallel.pool.Constant(solver_type);

%% Wavy wall features------------------------------------------------------
%Ny = 48; Nx = 36; Nz = 4; % Nz to Nx for streamwise terms
Ny = 120; Nx = 96; Nz = 4; % Nz to Nx for streamwise terms

epsilon = 0.50;
h = 1; y0 = h;
A1 = epsilon;
Lx = 9.0;
Ly = 1.0 + 1.0*epsilon;
Re = 100;

%% Plot periodic hill curve
x_val = linspace(0, Lx, Nx);
A  = 4.5;
B  = 3.5;
C  = 1/6;

y1_vals = 3*C * ( 1 + tanh( B*(abs(x_val - A) - B) ) );   % bottom wall y'(x,0)

% figure;
% plot(x_val,y1_vals)
% hold ons

%% Input from Dedalus
load data_x.mat
load data_y.mat
load mask_smooth_hillperiodic_120x96.mat
x_old = data_x; y_old = data_y;
%y_old = y_old + abs(data_y(1));

%% -------------fine mesh from Dedalus------------------------------------- 
[X_old, Y_old] = meshgrid(x_old, y_old);  
%% ------------------------------------------------------------------------

load u_mean_zt.mat
load v_mean_zt.mat
%load w_mean_zt.mat
load dUdx_mean_zt.mat
load dUdy_mean_zt.mat
load dVdx_mean_zt.mat
load dVdy_mean_zt.mat

U_f  = u_mean_zt;  V_f  = v_mean_zt;     % [Ny_f x Nx_f]
Ux_f = dUdx_mean_zt; Uy_f = dUdy_mean_zt; 
Vx_f = dVdx_mean_zt; Vy_f = dVdy_mean_zt; 
K_f  = mask_smooth;                        % [Ny_f x Nx_f]

% % interp2 requires increasing vectors; Chebyshev y often descending
% if numel(y_old)>1 && y_old(2) < y_old(1)
%     y_old = fliplr(y_old);
%     U_f  = flipud(U_f);  V_f  = flipud(V_f);
%     Ux_f = flipud(Ux_f); Uy_f = flipud(Uy_f);
%     Vx_f = flipud(Vx_f); Vy_f = flipud(Vy_f);
%     K_f  = flipud(K_f);
% end

%% Spectral operators
[x,Dxx] = fourdif(Nx,2);
[~,Dx]  = fourdif(Nx,1);
[yc,DM] = chebdif(Ny,2);  D2 = DM(:,:,2);  D1 = DM(:,:,1);

% Scaling
D1  = -D1*(2/(2+2*epsilon)); 
D2  = D2*(2/(2+2*epsilon))^2;
Dx  = Dx*(2*pi/Lx);
Dxx = Dxx*(2*pi/Lx)^2;

Ix = speye(Nx); Iy = speye(Ny);
I  = kron(Ix,Iy);  % sparse because Ix, Iy are sparse

[~,w] = clencurt(Ny-1);

% Scaling
x = x * (Lx / (2*pi));
y = yc/2*(2+2*epsilon);
y = flip(y); % chebyshev is in descending order so needs to be flipped
w = w*((2+2*epsilon)/2); % normalize weighted matrix

% w = spdiags(reshape(w,[],1),0,Ny,Ny);
% w = kron(Ix, diag(w));     % weighted matrix
w_sqrt = sqrt(w); 
w_sqrt = spdiags(w_sqrt',0,Ny,Ny);
w_sqrt = kron(Ix, w_sqrt);

w_sqrt_inv = 1./(sqrt(w)); % elementwise sqrt
w_sqrt_inv = spdiags(w_sqrt_inv',0,Ny,Ny);
w_sqrt_inv = kron(Ix, w_sqrt_inv);

x_old = reshape(x_old, 1, []);   % force row vector
y_old = reshape(y_old, [], 1);   % force column vector

%% Coarse grid
[X_new, Y_new] = meshgrid(x, y);

%% (1) Interpolate (vector form) to coarse grid
U_xy_coarse   = interp2(X_old, Y_old, U_f,  X_new, Y_new, 'linear');
V_xy_coarse   = interp2(X_old, Y_old, V_f,  X_new, Y_new, 'linear');
dUdx_xy_coarse= interp2(X_old, Y_old, Ux_f, X_new, Y_new, 'linear');
dUdy_xy_coarse= interp2(X_old, Y_old, Uy_f, X_new, Y_new, 'linear');
dVdx_xy_coarse= interp2(X_old, Y_old, Vx_f, X_new, Y_new, 'linear');
dVdy_xy_coarse= interp2(X_old, Y_old, Vy_f, X_new, Y_new, 'linear');

% % U_xy_coarse   = interp2(x_old, y_old, U_f,  X_new, Y_new, 'linear');
% % V_xy_coarse   = interp2(x_old, y_old, V_f,  X_new, Y_new, 'linear');
% % dUdx_xy_coarse= interp2(x_old, y_old, Ux_f, X_new, Y_new, 'linear');
% % dUdy_xy_coarse= interp2(x_old, y_old, Uy_f, X_new, Y_new, 'linear');
% % dVdx_xy_coarse= interp2(x_old, y_old, Vx_f, X_new, Y_new, 'linear');
% % dVdy_xy_coarse= interp2(x_old, y_old, Vy_f, X_new, Y_new, 'linear');

%% Flattened fields (keep your names)
N = Ny*Nx;
c_number = 12;

U_vec = reshape(U_xy_coarse,[],1);
V_vec = reshape(V_xy_coarse,[],1);

U = spdiags(U_vec,0,N,N);              % (2) sparse diag
V = spdiags(V_vec,0,N,N);

dUdx = spdiags(reshape(dUdx_xy_coarse,[],1),0,N,N);
dUdy = spdiags(reshape(dUdy_xy_coarse,[],1),0,N,N);
dVdx = spdiags(reshape(dVdx_xy_coarse,[],1),0,N,N);
dVdy = spdiags(reshape(dVdy_xy_coarse,[],1),0,N,N);

x_new = kron(diag(x),Iy);              % keep for any plotting (small)
y_new = kron(diag(y),Ix);

% Penalty mask
K = spdiags(reshape(mask_smooth,[],1),0,N,N);  % (2) sparse diag

%% Wavenumbers/frequency
kxn = 1; kzn = 36; kx = 1;
omega = 1;
kx_list = logspace(-4,0.48,kxn);
kz_list = logspace(-2,1.2,kzn);
c_list  = linspace(-1,1,c_number);

%% Allocate once (2) sparse
A = spalloc(4*N,4*N, 50*4*N);                 % rough nnz guess
B = spalloc(4*N,3*N, 3*N);
C = spalloc(3*N,4*N, 3*N);
E = spalloc(4*N,4*N, 3*N);
D = spalloc(3*N,3*N, 3*N);
D_new = spalloc(9*N,3*N, 9*N);
A_init_new = spalloc(4*N,4*N, 4*N);

zero_matrix = spalloc(N,N, N);
B_tilde = [w_sqrt_inv zero_matrix zero_matrix; 
           zero_matrix w_sqrt_inv zero_matrix;
           zero_matrix zero_matrix w_sqrt_inv;
           zero_matrix zero_matrix zero_matrix];

C_tilde = [w_sqrt zero_matrix zero_matrix zero_matrix;
           zero_matrix w_sqrt zero_matrix zero_matrix;
           zero_matrix zero_matrix w_sqrt zero_matrix];

E = [I zero_matrix zero_matrix zero_matrix;
     zero_matrix I zero_matrix zero_matrix;
     zero_matrix zero_matrix I zero_matrix;
     zero_matrix zero_matrix zero_matrix zero_matrix];

% Boundary row indices (3) vectorized
topRows    = (0:Nx-1)*Ny + 1;
bottomRows = (1:Nx)*Ny;
bd = [topRows, bottomRows];

  B_tilde(bd, :)        = 0;  B_tilde(bd+N, :) = 0;  B_tilde(bd+2*N, :) = 0;
  E(bd, :)        = 0;  E(bd+N, :) = 0;  E(bd+2*N, :) = 0;


% ---- Robust parallel setup: prefer threads; safe up to 4 workers ----
delete(gcp('nocreate'));

% Respect SLURM cpus-per-task; cap at 4
slurm_cpt = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(slurm_cpt) || isempty(slurm_cpt), slurm_cpt = 8; end
nworkers = max(1, min(4, slurm_cpt));

% Prevent nested threading inside each worker (important!)
setenv('OMP_NUM_THREADS', '1');
setenv('MKL_NUM_THREADS', '1');
maxNumCompThreads(nworkers);   % optional: hint MATLAB BLAS limit

try
    parpool('threads', nworkers);
catch
    warning('Thread pool unavailable; continuing in serial.');
end

tic
for i = 20  % pick a kz index (unchanged logic)
    kz = kz_list(i);

    % (5) cache per-i pieces
    laplacian = -kz^2*kron(Ix,Iy) + kron(Ix,D2) + kron(Dxx,Iy); % sparse

    % Interior blocks (sparse)
    A11 = U*kron(Dx,Iy) + V*kron(Ix,D1) + dUdx - (1/Re).*laplacian;
    A12 = dUdy;
    A21 = dVdx;
    A14 = kron(Dx,Iy);
    A22 = U*kron(Dx,Iy) + V*kron(Ix,D1) + dVdy - (1/Re)*laplacian;
    A24 = kron(Ix,D1);
    A33 = U*kron(Dx,Iy) + V*kron(Ix,D1) - (1/Re)*laplacian;
    A34 = 1i*kz*I;
    A41 = kron(Dx,Iy);
    A42 = kron(Ix,D1);
    A43 = 1i*kz*kron(Ix,Iy);

    % Assemble A (sparse)
    A = -[A11 A12 sparse(N,N) A14; ...
          A21 A22 sparse(N,N) A24; ...
          sparse(N,N) sparse(N,N) A33 A34; ...
          A41 A42 A43 sparse(N,N)];

    % B, C, E (sparse identities in their blocks) (2)
%     B = spalloc(4*N,3*N, 3*N);
%     B(1:N,1:N)             = I;
%     B(N+1:2*N,N+1:2*N)     = I;
%     B(2*N+1:3*N,2*N+1:3*N) = I;
% 
%     E = spalloc(4*N,4*N, 3*N);
%     E(1:N,1:N)             = I;
%     E(N+1:2*N,N+1:2*N)     = I;
%     E(2*N+1:3*N,2*N+1:3*N) = I;
% 
%     C = spalloc(3*N,4*N, 3*N);
%     C(1:N,1:N)             = I;
%     C(N+1:2*N,N+1:2*N)     = I;
%     C(2*N+1:3*N,2*N+1:3*N) = I;

    % (3) Vectorized boundary enforcement on A, and zero same rows in B & E
    bd_all = [bd, bd+N, bd+2*N];  % u,v,w rows in 4N system
    % zero rows then set diagonal entry to 1 on those rows & their variable blocks
    A(bd, :)        = 0;  A(sub2ind([4*N,4*N], bd, bd)) = 1;
    A(bd+N, :)      = 0;  A(sub2ind([4*N,4*N], bd+N, bd+N)) = 1;
    A(bd+2*N, :)    = 0;  A(sub2ind([4*N,4*N], bd+2*N, bd+2*N)) = 1;

%     B(bd, :)        = 0;  B(bd+N, :) = 0;  B(bd+2*N, :) = 0;
%     E(bd, :)        = 0;  E(bd+N, :) = 0;  E(bd+2*N, :) = 0;

    % Penalty (2) sparse blkdiag applied once
    K_block = blkdiag(K, K, K, sparse(N,N));
    A = A - K_block;

%     % Weighting / tilde blocks (unchanged)
%     C_tilde = blkdiag(w.^0.5, w.^0.5, w.^0.5) * C;
%     D = blkdiag(w.^0.5, w.^0.5, w.^0.5);
%     B_tilde = B / D;

%     % keep sparse
%     A = sparse(A);
%     B_tilde = sparse(B_tilde);
%     C_tilde = sparse(C_tilde);
%     E = sparse(E);

    % cache per i (5)
    A_s = A; E_s = E; B_s = B_tilde; C_s = C_tilde;
    solver = lower(solver_type);

    omega = zeros(c_number,1);
    result_sigma = zeros(c_number, kzn);
    result_U_svd = cell(c_number, kzn);
    result_V_svd = cell(c_number, kzn);

    opts.tol   = 1e-8;
    opts.maxit = 200;
    opts.p     = 5;

    tStart = tic;
    eigen_size = 1;

    % Share large, read-only matrices to avoid copies in workers
    A_c = parallel.pool.Constant(A_s);
    E_c = parallel.pool.Constant(E_s);
    B_c = parallel.pool.Constant(B_s);
    C_c = parallel.pool.Constant(C_s);

    parfor c_index = 1:c_number
    sigma1 = NaN; Hv1 = [];

    % pull shared matrices from Constants (no worker copies)
    A_loc = A_c.Value;  E_loc = E_c.Value;  B_loc = B_c.Value;  C_loc = C_c.Value;

    omega(c_index) = -c_list(c_index) * kx;
    omega_val = omega(c_index);

    switch solver
        case 'svds'
            M   = decomposition(1i*omega_val*E_loc - A_loc, 'lu');
            M_T = decomposition((1i*omega_val*E_loc - A_loc)', 'lu');

            Aop    = @(x) M \ (B_loc * x);
            Aadjop = @(x) (B_loc' * (M_T \ (C_loc' * x)));
            Hhandle = @(x,flag) Hfun_wrapper_LU(x, flag, Aop, Aadjop, C_loc);

            m = size(C_loc,1);
            n = size(B_loc,2);

            [U_svd, S2_svd, V_svd] = svds(Hhandle, [m,n], 1, 'largest', opts);
            U1_svd = U_svd;  V1_svd = V_svd;  sigma1 = S2_svd(1,1);

            Hv1 = C_loc * Aop(V1_svd);

        case 'svd'
            H = C_loc * ((1i*omega_val*E_loc - A_loc) \ B_loc);
            [U_svd, S2_svd, V_svd] = svd(H, 'econ');
            U1_svd = U_svd(:,1); V1_svd = V_svd(:,1); sigma1 = S2_svd(1,1);

            Hv1 = C_loc * ((1i*omega_val*E_loc - A_loc) \ (B_loc * V1_svd));
    end

    rel_residual = norm(Hv1 - sigma1 * U1_svd) / max(norm(Hv1), eps);
    fprintf("omega = %.4f, Sigma = %0.4e, Residual = %.2e\n", omega_val, sigma1, rel_residual);

    result_sigma(c_index, i) = sigma1;
    result_U_svd{c_index, i} = U1_svd;
    result_V_svd{c_index, i} = V1_svd;

    % --- free large temporaries to lower peak memory on workers ---
    M   = []; M_T = []; U_svd = []; V_svd = []; S2_svd = []; Hv1 = [];
end


    elapsed = toc(tStart);
    fprintf("Elapsed time for i = %d is %.2f seconds\n", i, elapsed);
    i
end

U_hat  = zeros(Ny, Nx, c_number, 'single');
V_hat  = zeros(Ny, Nx, c_number, 'single');
W_hat  = zeros(Ny, Nx, c_number, 'single');

U2_hat = zeros(Ny, Nx, c_number, 'single');
V2_hat = zeros(Ny, Nx, c_number, 'single');
W2_hat = zeros(Ny, Nx, c_number, 'single');

for c_index = 1:c_number
    vec = result_U_svd{c_index, i};  assert(numel(vec) == 3*N); % divide by weighting D matrix
    vec2= result_V_svd{c_index, i};  assert(numel(vec2)== 3*N);

    % RESPONSE
    U_hat(:,:,c_index)  = reshape(w_sqrt_inv*real(vec(1:N)),            Ny, Nx);
    V_hat(:,:,c_index)  = reshape(w_sqrt_inv*real(vec(N+1:2*N)),        Ny, Nx);
    W_hat(:,:,c_index)  = reshape(w_sqrt_inv*real(vec(2*N+1:3*N)),      Ny, Nx);

    % FORCING
    U2_hat(:,:,c_index) = reshape(w_sqrt_inv*real(vec2(1:N)),           Ny, Nx);
    V2_hat(:,:,c_index) = reshape(w_sqrt_inv*real(vec2(N+1:2*N)),       Ny, Nx);
    W2_hat(:,:,c_index) = reshape(w_sqrt_inv*real(vec2(2*N+1:3*N)),     Ny, Nx);
end
elapsed = toc;
disp(elapsed);

save('stability_results.mat', ...
     'kz_list', 'x', 'y', ...
     'omega', 'c_list', ...
     'U_hat', 'V_hat', 'W_hat', ...
     'U2_hat', 'V2_hat', 'W2_hat','result_sigma', "result_U_svd", "result_V_svd");
