%% main signed distance function , binary mask and normalized mask for input putput SVDS
clc; clear all; close all
% Define domain
%Lx = 0.6 * pi;

% %% Read dedalus data file for x y dataset for grid and velocity profile
% read_dedalus_data;
% %% ---------------------xxxxxxxxxxxxxxxxxxxxxxxxx--------------------------
%folderpath = 'E:/dedalus/Wavy Wall/dedalus_local_Re75_3/snapshots_channel/mean_v_Re75.00_c1.00' 
%%folderpath = 'E:\dedalus\Wavy Wall\flexible signed distance function\opt_erf_peclet\dedalus_local_wavy_wall_Re290_a1\snapshots_channel\mean_v_Re290.00_c1.00';
folderpath = 'E:/dedalus/Wavy Wall/flexible signed distance function/opt_erf_peclet/dedalus_local_Re190_nek5000_hill_ubulk_2D/snapshots_channel/mean_v_Re190.00_c1.00'

cd(folderpath)
%Lx = 4.3;
%epsilon = 0.54; %0.125
%Ly = 1.0 + epsilon;

load data_x.mat
load data_y.mat

%Ny = 120; Nx = 96;  % Nz to Nx for streamwise terms
Ny = 140; Nx = 108;  % Nz to Nx for streamwise terms

% Wavy validation case Choo parameters
epsilon = 0.50;
%epsilon = 0.25;
h = 1;
y0 = h;

A1 = epsilon;
A2 = epsilon;
%Lx = 0.6 * pi;
Lx = 9.0;
Ly = 1.0 + 1.1*epsilon;

%Re = 75.30;
Re = 190;

[x,Dx] = fourdif(Nx,1); % first derivative x
[y_cheb,DM] = chebdif(Ny,2);
y = y_cheb*Ly;
y = flip(y);
x = x * (Lx / (2*pi)); 

[x_grid, y_grid] = meshgrid(x, y);

% Wall parameters
dy = Ly / Ny;
%mask_threshold = 2 * 0.5 * dy; % 8 to 6 to 3
mask_option = 6;

c = 1.0; 
mask_const = Re/c;  % Example value for penalization stiffness: 400
eta       = c*(1.0/Re)
gamma = sqrt(eta / Re);                      % ε = sqrt(η / Re)
mask_threshold = 3.11346786 * gamma;         % δ* = 3.113... * ε

%% Compute distance and mask
%[d_perp] = compute_signed_distance_mask(x_grid, y_grid, Lx, y0, A1, A2);
[d_perp] = compute_signed_distance_mask_nek5000(x_grid, y_grid, Lx, -y0, A1);

% Construct binary mask: solid = 1, interface = 0.5, fluid = 0
    binary_mask = zeros(size(d_perp));
    solid_region = d_perp < -mask_threshold;
    fluid_region = d_perp >  mask_threshold;
    interface_region = ~(solid_region | fluid_region);

    binary_mask(solid_region) = 1;
    binary_mask(interface_region) = 0.5;
  
if mask_option == 1
    %% Normalized smooth mask using erf
    [mask_smooth,mask_smooth_solid] = normalized_mask_erf(d_perp, mask_threshold, mask_const);

elseif mask_option == 2
    %% Normalized and shifted smooth mask using erf
    [mask_smooth, mask_smooth_solid, mask_smooth_fluid] = normalized_mask1_shifted(d_perp, mask_threshold, mask_const);

elseif mask_option == 3
    %% Normalized optimized smooth mask using erf
    [mask_smooth, mask_smooth_solid, mask_smooth_fluid, mask_smooth_interface] = normalized_mask_optimized(d_perp,Re, c, mask_const);

elseif mask_option == 4
     %% Normalized optimized smooth mask using tanh
    [mask_smooth, mask_smooth_solid, mask_smooth_fluid] = normalized_mask_optimized_tanh(d_perp, mask_const, Re)
elseif mask_option == 5
     %% Normalized optimized smooth mask using optimized erf
    [mask_smooth, mask_smooth_solid, mask_smooth_fluid, mask_smooth_interface] = normalized_mask_optimized3(d_perp, Re, mask_const, c, Ny, Nx, Lx, Ly)

elseif mask_option == 6
     %% Normalized optimized smooth mask using peclet number and optimized erf
     Hav = y0;
     Pe_max = 2.0;
     kappa = 1.0;
    [mask_smooth, mask_smooth_combined, mask_smooth_solid, mask_smooth_fluid, mask_smooth_interface, delta_star, delta_floor, delta_opt, delta_burns, delta_pe, Pe_delta, s] = ...
    normalized_mask_optimized_peclet(d_perp, Re, mask_const, c, Ny, Ly, Lx, Nx, Hav, Pe_max, kappa,data_x,data_y)

end

%% Normalized mask using tanh function
mask_smooth_tanh = normalized_mask_tanh(d_perp, mask_threshold, mask_const);

%% Coordinates of solid region from mask_smooth_solid
[x_solid, y_solid] = deal(x_grid(mask_smooth_solid), y_grid(mask_smooth_solid));
[x_interface, y_interface] = deal(x_grid(mask_smooth_interface), y_grid(mask_smooth_interface));

solid_nodes = [x_solid y_solid];

%% Coordinates of solid region from mask_smooth_solid
[x_fluid, y_fluid] = deal(x_grid(mask_smooth_fluid), y_grid(mask_smooth_fluid));

fluid_nodes = [x_fluid y_fluid];

%% Find no slip region in interface
% x_wall2 = data_x;
% y_wall2 = data_y;
% y_wall2a = -y0 - A1 * sin(2 * pi / Lx * data_x);
% y_wall2b = y0 + A1 * sin(2 * pi / Lx * data_x);


% % --- dense analytic wall (352 pts = Ny) ---
% x_wall2   = linspace(0, Lx, Ny);           % 352 samples in x
% y_wall2a_dense = -y0 - A1 .* sin(2*pi/Lx .* x_wall2);
% y_wall2b_dense =  y0 + A2 .* sin(2*pi/Lx .* x_wall2);
% 
% y_wall2a = y_wall2a_dense;
% y_wall2b = y_wall2b_dense;

% --- values at your simulation x grid (Nx=288) ---
% y_wall2a = interp1(x_wall2, y_wall2a_dense, data_x, 'pchip');  % or 'linear'
% y_wall2b = interp1(x_wall2, y_wall2b_dense, data_x, 'pchip');

% Optional: visualize the wall
% figure;
% scatter(x_wall, y_wall, 5, 'r', 'filled');
% xlabel('x'); ylabel('y');
% title('Wall Surface: Points where d\_perp ≈ 0');
% axis equal;

%% Post Processing
% Plot1
% figure(1);
% pcolor(x_grid, y_grid, binary_mask); shading flat; axis equal tight;
% colorbar; title('Binary Mask (Solid=1, Interface=0.5, Fluid=0)');

% Plot2
figure(2);
pcolor(x_grid, y_grid, mask_smooth); shading flat; axis equal tight;
colorbar; title(sprintf('Smooth Normalized Mask (erf) at Re = %.2f', Re));
outdir      = folderpath;  % or wherever you like
plot_fname  = sprintf('smooth_mask_Re%.2f_mask%d.png', Re, mask_const);
full_plotfn = fullfile(outdir, plot_fname);
saveas(gcf, full_plotfn);

% plot_filename = sprintf('smooth normalied mask(erf) for Re%.2f and mask %d.png', Re, mask_const);
% saveas(gcf, plot_filename);

% % Plot3
% figure(3);
% pcolor(x_grid, y_grid, mask_smooth_tanh); shading flat; axis equal tight;
% colorbar; title('Smooth Normalized Mask (tanh)');

% Plot4
figure(4)
plot(x_solid,y_solid,'*');
hold on
% Plot5
plot(x_fluid,y_fluid,'o');
title("solid and fluid cells")
legend('solid cells','fluid cells','interface')
xlim([0, Lx]);
ylim([-Ly, Ly]);
hold on

  %dlmwrite('mask_smooth_erf.csv',mask_smooth); 
% dlmwrite('mask_smooth_tanh.csv',mask_smooth_tanh);
% dlmwrite('solid_nodes.csv',solid_nodes); 

% Output files
  %save("mask_smooth_hillperiodic_Re100_132x108.mat", "mask_smooth");
 % save("mask_smooth_hpc_80x64.mat");

filename = ['mask_smooth_hillperiodic_Re' num2str(Re) '_' num2str(Ny) 'x' num2str(Nx) '.mat'];
save(filename, 'mask_smooth');


  % save("mask_smooth_fluid");
  save("mask_smooth_solid.mat");
  save("x_solid.mat","x_solid");
  save("y_solid.mat","y_solid");


function [mask_smooth, mask_smooth_solid] = normalized_mask_erf(d_perp, mask_threshold, mask_const)
    % Compute a smooth penalization mask based on signed distance
    % - Solid: mask_const
    % - Fluid: 0
    % - Interface: erf transition

    steepness = 1;  % Can increase for sharper interface (e.g., 2), or decrease for smoother

    mask_smooth = zeros(size(d_perp));

    % Region masks
    solid_region     = (d_perp <= -mask_threshold);
    fluid_region     = (d_perp >  mask_threshold);
    interface_region = ~solid_region & ~fluid_region;
        
    mask_smooth = 0.5 * (1 - erf(steepness * sqrt(pi) * d_perp(interface_region) / mask_threshold)) * mask_const;
    % Output solid region mask (binary, same shape)
    mask_smooth_solid = solid_region;
    
end


function mask_smooth = normalized_mask_tanh(d_perp, mask_threshold, mask_const)
    % Compute smooth mask using tanh transition (alternative to erf)
    % - Solid: mask_const
    % - Fluid: 0
    % - Interface: smooth transition via tanh

    mask_smooth = zeros(size(d_perp));

    % Region masks
    solid_region     = (d_perp <= -mask_threshold);
    fluid_region     = (d_perp >  mask_threshold);
    interface_region = ~solid_region & ~fluid_region;

    % Assign values
    mask_smooth(solid_region) = mask_const;
    mask_smooth(fluid_region) = 0;
    mask_smooth(interface_region) = ...
        0.5 * (1 - tanh(d_perp(interface_region) / mask_threshold)) * mask_const;
end

function [mask_smooth, mask_smooth_solid, mask_smooth_fluid] = normalized_mask1_shifted(d_perp, mask_threshold, mask_const)
%NORMALIZED_MASK1_SHIFTED   Shifted erf mask + region indicators (O(ε²) error)
%
%   [mask_smooth, mask_smooth_solid, mask_smooth_fluid] = ...
%       normalized_mask1_shifted(d_perp, mask_threshold, mask_const)
%
%   Inputs:
%     d_perp         Array of signed distances from the wall
%     mask_threshold Half‑width δ of the original erf ramp
%     mask_const     Penalty strength (1/τ)
%
%   Outputs:
%     mask_smooth        Shifted erf mask array
%     mask_smooth_solid  Logical array: true inside the solid region
%     mask_smooth_fluid  Logical array: true inside the fluid region

    % 1) compute displacement length ℓ = δ / √π
    shift = mask_threshold / sqrt(pi);

    % 2) shift distance deeper into the solid
    d_s = d_perp + shift;

    % 3) define regions on the shifted coordinate
    solid_region = (d_s <= -mask_threshold);
    fluid_region = (d_s >=  mask_threshold);
    interface_region  = ~(solid_region | fluid_region);

    % 4) build the shifted erf mask
    % 4) smooth mask (not piecewise anymore)
    mask_smooth = 0.5 * (1 - erf(sqrt(pi) * d_s / mask_threshold)) * mask_const;

    mask_smooth_solid = solid_region;
    mask_smooth_fluid = fluid_region;
    %mask_smooth_interface = interface_region;

end
