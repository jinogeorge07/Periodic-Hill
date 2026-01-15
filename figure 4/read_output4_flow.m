clc; clear all;close all;

folderpath = './'
%=========================================================================================================================

mask_axis = "False"
data_x = readmatrix('data_x.csv');
data_y = readmatrix('data_y.csv');   

u_xy = readmatrix("u_xy.csv")
mask_smooth = readmatrix('mask_smooth.csv');   
snap_idx = readmatrix('snap_idx.csv'); 

Re = 100;
Nt = 101;
% For RE = 190
if Re == 190
    Ny = 576; % 48, 80
    Nx = 384; % 36, 64
else
    % For RE = 100
    Ny = 512; % 48, 80
    Nx = 352; % 36, 64
    Nz = 2;
end

% v1_flat = readmatrix('v1.csv');   
% v2_flat = readmatrix('v2.csv');   
% 
% % v1 = reshape(v1_flat, Nx, Ny, Nt);
% % v2 = reshape(v2_flat,Nx,Ny,Nt);
% 
% v1 = reshape(v1_flat, Ny, Nx, Nt);
% v2 = reshape(v2_flat, Ny, Nx, Nt);
% 
% v1 = permute(reshape(v1_flat(:), Ny, Nx, Nt), [2 1 3]);
% v2 = permute(reshape(v2_flat(:), Ny, Nx, Nt), [2 1 3]);

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

x_val = linspace(0, Lx, Nx);  % replaces fourdif x values for plotting only
y1_vals =  -y0 + A1*tanh(3.5*abs(x_val-4.5)- 3.5)
y_wall2a = -y0 + A1 *tanh(3.5*abs(data_x-4.5)- 3.5);


% Pick centerline slices
n = Nx;
x_idx_vec = linspace(1,n,n);  % <-- EDIT THIS VECTOR
z_idx = 1;          % center in z
t_idx = snap_idx(end);        % last snapshot
u_thresh = 1e-3;              % velocity threshold
x_tol = 1e-6;                 % tolerance for wall match
%mask_smooth = squeeze(mask_smooth(:,:,t_idx));

%v1 = permute(v1, [2, 1, 3]);  % reorder to (x, y, z, t)
%v2 = permute(v2, [2, 1, 3]);  % reorder to (x, y, z, t)

%mask_smooth = mask_smooth';     % transpose if needed
l_star_grid = nan(length(x_idx_vec), 2);  % preallocate (2 l* per x)

%mask_midline_hillp;

%ywall_values;
mask_smooth = permute(mask_smooth,[2,1]);
u_xy = permute(u_xy,[2,1]);

% Create output directory for plots
output_dir = sprintf('plots_Re_%g_mask_%.1f', Re);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


for k = 1:length(x_idx_vec)
    x_idx = x_idx_vec(k);

    % Extract velocity and mask
%     u_y = squeeze(v1(x_idx, :, t_idx));
%     v_y = squeeze(v2(x_idx, :, t_idx));
    u_y = u_xy(x_idx,:);
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
    filename_safe = sprintf('u velocity_mask_x_%.3f_t_%d_mask_%.1f.png', data_x(x_idx), t_idx, Re);
    saveas(gcf, fullfile(output_dir, filename_safe));


    %% end of plot --------------xxxxxxxxxxxxxxxxxx------------------------
    end


end