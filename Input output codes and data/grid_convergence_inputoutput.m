clc; close all;clear all;


folderpath = 'E:/dedalus/Wavy Wall/flexible signed distance function/opt_erf_peclet/dedalus_local_Re100_nek5000_hill_2a_ubulk_2D/snapshots_channel/mean_v_Re100.00_c1.00/grid convergence_kz5'
cd(folderpath)

%% Wavenumbers/frequency
c_number = 12;
kxn = 1; kzn = 36; kx = 1;
omega = 1;
kx_list = logspace(-4,0.48,kxn);
kz_list = logspace(-2,1.2,kzn);
c_list  = linspace(-1,1,c_number);
omega = -c_list * kx;


load max_Uhat_60x48.mat
load max_Vhat_60x48.mat

load max_Uhat_72x60.mat
load max_Vhat_72x60.mat

load max_Uhat_96x80.mat
load max_Vhat_96x80.mat

load max_Uhat_120x96.mat
load max_Vhat_120x96.mat

load max_Uhat_132x108.mat
load max_Vhat_132x108.mat

% u_diff1 = abs(max_Uhat_120x96 - max_Uhat_96x80)./max_Uhat_120x96;
% u_diff2 = abs(max_Uhat_120x96 - max_Uhat_72x60)./max_Uhat_120x96;
% u_diff3 = abs(max_Uhat_120x96 - max_Uhat_60x48)./max_Uhat_120x96;
%
% v_diff1 = abs(max_Vhat_120x96 - max_Vhat_96x80)./max_Vhat_120x96;
% v_diff2 = abs(max_Vhat_120x96 - max_Vhat_72x60)./max_Vhat_120x96;
% v_diff3 = abs(max_Vhat_120x96 - max_Vhat_60x48)./max_Vhat_120x96;


u_diff1 = abs(max_Uhat_132x108 - max_Uhat_120x96)./max_Uhat_132x108;
u_diff2 = abs(max_Uhat_132x108 - max_Uhat_96x80)./max_Uhat_132x108;
u_diff3 = abs(max_Uhat_132x108 - max_Uhat_72x60)./max_Uhat_132x108;
u_diff4 = abs(max_Uhat_132x108 - max_Uhat_60x48)./max_Uhat_132x108;


v_diff1 = abs(max_Vhat_132x108 - max_Vhat_120x96)./max_Uhat_132x108;
v_diff2 = abs(max_Vhat_132x108 - max_Vhat_96x80)./max_Vhat_132x108;
v_diff3 = abs(max_Vhat_132x108 - max_Vhat_72x60)./max_Vhat_132x108;
v_diff4 = abs(max_Vhat_132x108 - max_Vhat_60x48)./max_Vhat_132x108;

figure
plot(omega,u_diff1,'o-', 'LineWidth', 2, 'MarkerSize', 6)
hold on
plot(omega,u_diff2,'s--', 'LineWidth', 2, 'MarkerSize', 6)
hold on
plot(omega,u_diff3,'d-.', 'LineWidth', 2, 'MarkerSize', 6)
hold on
plot(omega,u_diff4,'o-', 'LineWidth', 2, 'MarkerSize', 6)
hold on
grid off;
%title('Resolvent Singular Values vs Phase Speed', 'FontSize', 22);
xlabel('Temporal frequency \omega', 'FontSize', 24);
ylabel('Relative error', 'FontSize', 24);

legend({ ...
    sprintf('(Grid5 - Grid4)/Grid5'), ...
    sprintf('(Grid4 - Grid3)/Grid5'), ...
    sprintf('(Grid3 - Grid2)/Grid5'), ...
    sprintf('(Grid2 - Grid1)/Grid5')}, ...
    'FontSize', 18, 'Location', 'northwest');

% Set axis tick label font size
set(gca, 'FontSize', 24);
lgd = legend;
set(lgd, 'Color', 'none');   % remove background
set(lgd, 'Box', 'off');      % optional: remove border
xlim([-1 1])
ylim([0 1])

saveas(gcf, 'resolvent_spectrum_deviation_U.png');


figure
plot(omega,v_diff1,'o-', 'LineWidth', 2, 'MarkerSize', 6)
hold on
plot(omega,v_diff2,'s--', 'LineWidth', 2, 'MarkerSize', 6)
hold on
plot(omega,v_diff3,'d-.', 'LineWidth', 2, 'MarkerSize', 6)
hold on
plot(omega,v_diff4,'o-', 'LineWidth', 2, 'MarkerSize', 6)
hold on
grid off;
%title('Resolvent Singular Values vs Phase Speed', 'FontSize', 22);
xlabel('Temporal frequency \omega', 'FontSize', 20);
ylabel('Relative error', 'FontSize', 20);

legend({ ...
    sprintf('(Grid1 - Grid2)/Grid1'), ...
    sprintf('(Grid1 - Grid3)/Grid1'), ...
    sprintf('(Grid1 - Grid4)/Grid1'), ...
    sprintf('(Grid1 - Grid5)/Grid1')}, ...
    'FontSize', 18, 'Location', 'northwest');

% Set axis tick label font size
set(gca, 'FontSize', 24);
lgd = legend;
set(lgd, 'Color', 'none');   % remove background
set(lgd, 'Box', 'off');      % optional: remove border
xlim([-1 1])
ylim([0 1])

saveas(gcf, 'resolvent_spectrum_deviation_V.png');



minimum_grid = [min(u_diff4) min(u_diff3) min(u_diff2) min(u_diff1)];
grid = [1 2 3 4]

figure();
plot(grid, minimum_grid,'s--', 'LineWidth', 2, 'MarkerSize', 6)
hold on
xlabel('Grid (coarse \rightarrow refined)', 'FontSize', 16);
ylabel('Relative error', 'FontSize', 16);
xlim([1 4])
ylim([0 0.5])
ax = gca;
ax.XTick = 1:4;
ax.XTickLabel = {'1','2','3','4'};
% Increase tick font size
ax.FontSize = 24;    % adjust as desired (16–22 is typical)
ax.TickLength = [0.02 0.02];   % optional: longer ticks
ax.LineWidth = 1.5;            % optional: thicker axes

hold on

saveas(gcf, 'grid convergence.png');
