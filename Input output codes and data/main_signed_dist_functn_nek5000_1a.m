%% main signed distance function , binary mask and normalized mask
%%clc; clear all; close all
% Define domain
%Lx = 0.6 * pi;

read_data = "false"
if read_data == "true"
    %% Read dedalus data file for x y dataset for grid and velocity profile
    read_dedalus_data;
    %% ---------------------xxxxxxxxxxxxxxxxxxxxxxxxx--------------------------
end

% Ny = 304; % 60, 80, 96
% Nx = 256; % 48, 64,64

% [x,Dx] = fourdif(Nx,1); % first derivative x
% [y_cheb,DM] = chebdif(Ny,2);
% y = y_cheb*Ly;
%
% x = x * (Lx / (2*pi));
x = data_x;
y = data_y;

[x_grid, y_grid] = meshgrid(x, y);

% Wall parameters
%y0 = 0.6;
y0 = 1.0;
A1 = epsilon;
A2 = epsilon;
%Re = 75.30
% y0 = 1.0;
% A1 = 0.125;
% A2 = 0.125;
%Re = 63.36;
dy = Ly / Ny;
%mask_threshold = 2 * 0.5 * dy; % 8 to 6 to 3
mask_option = 7;
distance_fn = "flexible"
Re
c = 1; alpha = 1.0
mask_const = Re^alpha/c;  % Example value for penalization stiffness: 400

if distance_fn == "const"
    mask_threshold = 6*0.5*dy

elseif distance_fn == "optimum"
    eta       = c*(1.0/Re)
    gamma = sqrt(eta / Re);                      % ε = sqrt(η / Re)
    delta_star = 3.11346786 * gamma;         % δ* = 3.113... * ε
    mask_threshold = 2*delta_star;
    %%mask_threshold = 2.0 * (3.11346786*0.01) % testing

elseif distance_fn == "flexible"
    % y is a 1-D vector of y-coordinates
    y = y(:);  % ensure column

    if numel(y) >= 2
        dy_row = zeros(size(y));
        dy_row(2:end-1) = 0.5*(y(3:end) - y(1:end-2));  % centered interior
        dy_row(1)       =  y(2) - y(1);                 % forward at bottom
        dy_row(end)     =  y(end) - y(end-1);           % backward at top
        mask_threshold  = 2.0 * max(abs(dy_row));       % scalar
    else
        mask_threshold  = 0;  % fallback if fewer than 2 points
    end

    % Chebyshev wall spacing
    %     j = 0:(Ny-1);
    %     y_cheb  = cos(pi * j / (Ny - 1));
    %     Y2_cheb = (Ly/2) .* y_cheb;
    %
    %     dy_wall = abs(Y2_cheb(2) - Y2_cheb(1));
    %     dx_wall = Lx / Nx;
    %     h_wall  = min(dx_wall, dy_wall);
    %     mask_threshold = h_wall;
    %     fprintf('dx_wall = %.16g\n', dx_wall);
    %     fprintf('dy_wall = %.16g\n', dy_wall);
    %     fprintf('h_wall = %.16g\n', h_wall);

end
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
    [mask_smooth,mask_smooth_solid, mask_smooth_fluid] = normalized_mask_erf(d_perp, mask_threshold, mask_const);

elseif mask_option == 2
    %% Normalized and shifted smooth mask using erf
    [mask_smooth, mask_smooth_solid, mask_smooth_fluid] = normalized_mask1_shifted(d_perp, mask_threshold, mask_const);

elseif mask_option == 3
    %% Normalized optimized smooth mask using erf
    [mask_smooth, mask_smooth_solid, mask_smooth_fluid, mask_smooth_interface] = normalized_mask_optimized(d_perp, Re, c, mask_const);

elseif mask_option == 4
    %% Normalized optimized smooth mask using erf
    [mask_smooth, mask_smooth_solid, mask_smooth_fluid, mask_smooth_interface] = normalized_mask_optimized_2(d_perp, Re, c, mask_const,Ny,Nx,Lx,Ly);

elseif mask_option == 5
    %% Normalized optimized smooth mask using tanh
    [mask_smooth, mask_smooth_solid, mask_smooth_fluid] = normalized_mask_optimized_tanh(d_perp, mask_const, Re)
elseif mask_option == 6
    %% Normalized optimized smooth mask using optimized erf
    [mask_smooth, mask_smooth_solid, mask_smooth_fluid, mask_smooth_interface] = normalized_mask_optimized3(d_perp, Re, mask_const, c, Ny, Nx, Lx, Ly)
elseif mask_option == 7
    %% Normalized optimized smooth mask using peclet number and optimized erf
    Hav = y0;
    alpha = 1.0;
    [mask_smooth, mask_smooth_combined, mask_smooth_solid, mask_smooth_fluid, mask_smooth_interface, delta_star, delta_floor, delta_opt, delta_burns, s] = ...
        normalized_mask_optimized_peclet(d_perp, Re, mask_const, c, Ny, Ly, Lx, Nx, Hav, alpha,data_x,data_y)

end

%% Normalized mask using tanh function
mask_smooth_tanh = normalized_mask_tanh(d_perp, mask_threshold, mask_const);

%% Coordinates of solid region from mask_smooth_solid
[x_solid, y_solid] = deal(x_grid(mask_smooth_solid), y_grid(mask_smooth_solid));

solid_nodes = [x_solid y_solid];

%% Coordinates of fluid region from mask_smooth_solid
%[x_fluid, y_fluid] = deal(x_grid(mask_smooth_combined), y_grid(mask_smooth_combined));
[x_fluid, y_fluid] = deal(x_grid(mask_smooth_fluid), y_grid(mask_smooth_fluid));

X_fluid = x_grid;
Y_fluid = y_grid;


% X_fluid(~mask_smooth_fluid) = NaN;
% Y_fluid(~mask_smooth_fluid) = NaN;

% if we need interface as well
X_fluid(~mask_smooth_combined) = NaN;
Y_fluid(~mask_smooth_combined) = NaN;

save('X_fluid.mat','X_fluid')
save('Y_fluid.mat','Y_fluid')

%% Find no slip region in interface

% % Define local grid size where we see max recirculation
% dy_local = abs(y(33) - y(32));
% tol = dy_local;  % or 1.5× if needed for robustness
% % Logical mask where signed distance ≈ 0 (i.e., wall surface)
% wall_region = abs(d_perp) <= tol;
%
% % Extract corresponding (x, y) coordinates from meshgrid
% x_wall = x_grid(wall_region);
% y_wall = y_grid(wall_region);

x_wall2 = data_x;
y_wall2 = data_y;
y_wall2a = -y0 + A1 *tanh(3.5*abs(data_x-4.5)- 3.5);
%y_wall2b = y0 + A1 * sin(2 * pi / Lx * data_x);


figure;
scatter(x_wall2, y_wall2a, 5, 'r', 'filled');
xlabel('x'); ylabel('y');
hold on

%% Post Processing
% Plot1
figure(1);
pcolor(x_grid, y_grid, binary_mask); shading flat; axis equal tight;
colorbar; title('Binary Mask (Solid=1, Interface=0.5, Fluid=0)');

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
xlim([0, 9]);
ylim([-1.5, 1.5]);
hold on

% dlmwrite('mask_smooth_erf.csv',mask_smooth);
% dlmwrite('mask_smooth_tanh.csv',mask_smooth_tanh);
% dlmwrite('solid_nodes.csv',solid_nodes);

% Output files
% save("mask_smooth.mat");
% save("mask_smooth_fluid");
save("mask_smooth_solid");

function [mask_smooth, mask_smooth_solid, mask_smooth_fluid] = normalized_mask_erf(d_perp, mask_threshold, mask_const)
% Compute a smooth penalization mask based on signed distance
% - Solid: mask_const
% - Fluid: 0
% - Interface: erf transition

steepness = 1;  % Can increase for sharper interface (e.g., 2), or decrease for smoother

%mask_smooth = zeros(size(d_perp));

% Region masks
solid_region     = (d_perp <= -mask_threshold);
fluid_region     = (d_perp >  mask_threshold);
interface_region = ~solid_region & ~fluid_region;

%mask_smooth = 0.5 * (1 - erf(steepness * sqrt(pi) * d_perp(interface_region) / mask_threshold)) * mask_const;
mask_smooth = 0.5 * (1 - erf(sqrt(pi) * d_perp / mask_threshold)) * mask_const;

% Output solid region mask (binary, same shape)
mask_smooth_solid = solid_region;
mask_smooth_fluid = fluid_region;
mask_smooth_interface = interface_region;

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
