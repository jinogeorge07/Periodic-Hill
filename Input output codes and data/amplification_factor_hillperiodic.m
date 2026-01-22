% Amplification factor and growth
clc; clear all; close all;
%folderpath = 'E:/dedalus/Wavy Wall/flexible signed distance function/opt_erf_peclet/dedalus_local_Re100_nek5000_hill_2a_ubulk_2D/snapshots_channel/mean_v_Re100.00_c1.00'
folderpath = 'E:/dedalus/Wavy Wall/flexible signed distance function/opt_erf_peclet/dedalus_local_Re190_nek5000_hill_ubulk_2D/snapshots_channel/mean_v_Re190.00_c1.00'
cd(folderpath)

%load('stability_results_Re100_120x96_kz20.mat');   % expects x,y,U_hat,V_hat,W_hat,kz_list,c_list,omega
Re = 190
c_number = 12;
compute_deviation = "false"
%% Wavenumbers/frequency
kxn = 1; kzn = 36; kx = 1;
omega = 1;
kx_list = logspace(-4,0.48,kxn);
kz_list = logspace(-2,1.2,kzn);
c_list  = linspace(-1,1,c_number);
omega = -c_list * kx;

filename = ['sigma_list_Re' num2str(Re) '_kz5.mat'];
data5 = load(filename);
filename = ['sigma_list_Re' num2str(Re) '_kz10.mat'];
data10 = load(filename);
filename = ['sigma_list_Re' num2str(Re) '_kz20.mat'];
data20 = load(filename);
filename = ['sigma_list_Re' num2str(Re) '_kz28.mat'];
data28 = load(filename);
filename = ['sigma_list_Re' num2str(Re) '_kz32.mat'];
data32 = load(filename);

sigma_list_5 = data5.sigma_list;
sigma_list_10 = data10.sigma_list;
sigma_list_20 = data20.sigma_list;
sigma_list_28 = data28.sigma_list;
sigma_list_32 = data32.sigma_list;

if compute_deviation == "true"

    data5a = load('sigma_list_Re50_kz5.mat');
    data10a = load('sigma_list_Re50_kz10.mat');
    data20a = load('sigma_list_Re50_kz20.mat');
    data32a = load('sigma_list_Re50_kz32.mat');

    sigma_list_5a = data5a.sigma_list;
    sigma_list_10a = data10a.sigma_list;
    sigma_list_20a = data20a.sigma_list;
    sigma_list_32a = data32a.sigma_list;


    amplification_kz5 = sigma_list_5- sigma_list_5a;
    amplification_kz10 = sigma_list_10- sigma_list_10a;
    amplification_kz20 = sigma_list_20- sigma_list_20a;
    amplification_kz32 = sigma_list_32- sigma_list_32a;

    plot(omega,amplification_kz5,'o-', 'LineWidth', 2, 'MarkerSize', 6)
    hold on
    plot(omega,amplification_kz10,'s--', 'LineWidth', 2, 'MarkerSize', 6)
    hold on
    plot(omega,amplification_kz20,'d-.', 'LineWidth', 2, 'MarkerSize', 6)
    hold on
    plot(omega,amplification_kz32,'o-', 'LineWidth', 2, 'MarkerSize', 6)
    hold on
    grid off;
    %title('Resolvent Singular Values vs Phase Speed', 'FontSize', 22);
    xlabel('Temporal frequency \omega', 'FontSize', 20);
    ylabel('Deviation in singular value \sigma', 'FontSize', 20);

    legend({ ...
        sprintf('\\Delta kz = %.2f', kz_list(5)), ...
        sprintf('\\Delta kz = %.2f', kz_list(10)), ...
        sprintf('\\Delta kz = %.2f', kz_list(20)), ...
        sprintf('\\Delta kz = %.2f', kz_list(32))}, ...
        'FontSize', 18, 'Location', 'northwest');

    % Set axis tick label font size
    set(gca, 'FontSize', 20);
    lgd = legend;
    set(lgd, 'Color', 'none');   % remove background
    set(lgd, 'Box', 'off');      % optional: remove border

    saveas(gcf, 'resolvent_spectrum_deviation.png');

end


sigma2_list_5  = sigma_list_5.^2;
sigma2_list_10 = sigma_list_10.^2;
sigma2_list_20 = sigma_list_20.^2;
sigma2_list_28 = sigma_list_28.^2;
sigma2_list_32 = sigma_list_32.^2;

% Amplification gain------------------------------------------------------
figure;
plot(c_list, sigma_list_5, 'o-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(c_list, sigma_list_10, 's--', 'LineWidth', 2, 'MarkerSize', 6);
plot(c_list, sigma_list_20, 'd-.', 'LineWidth', 2, 'MarkerSize', 6);
plot(c_list, sigma_list_28, 'd-.', 'LineWidth', 2, 'MarkerSize', 6);
plot(c_list, sigma_list_32, 'd-.', 'LineWidth', 2, 'MarkerSize', 6);

grid off;
%title('Resolvent Singular Values vs Phase Speed', 'FontSize', 22);
xlabel('Phase speed c', 'FontSize', 24);
ylabel('Singular value \sigma', 'FontSize', 24);

legend({ ...
    sprintf('kz = %.2f', kz_list(5)), ...
    sprintf('kz = %.2f', kz_list(10)), ...
    sprintf('kz = %.2f', kz_list(20)), ...
    sprintf('kz = %.2f', kz_list(28)), ...
    sprintf('kz = %.2f', kz_list(32))}, ...
    'FontSize', 18, ...
    'Orientation','horizontal', ...
    'NumColumns', 2, ...   % <-- first row: 3 items
    'Location','northwest');  % <-- inside, top-left

% Set axis tick label font size
set(gca, 'FontSize', 20);
lgd = legend;
set(lgd, 'Color', 'none');   % remove background
set(lgd, 'Box', 'off');      % optional: remove border
xlim([-1 1])
ylim([0 200])

%saveas(gcf, 'resolvent_spectrum_kz_all.png');
filename = ['resolvent_spectrum_kz_all_' num2str(Re) '.png'];
saveas(gcf, filename);

figure;
plot(omega, sigma_list_5, 'o-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(omega, sigma_list_10, 's--', 'LineWidth', 2, 'MarkerSize', 6);
plot(omega, sigma_list_20, 'd-.', 'LineWidth', 2, 'MarkerSize', 6);
plot(omega, sigma_list_28, 'd-.', 'LineWidth', 2, 'MarkerSize', 6);
plot(omega, sigma_list_32, 'p-.', 'LineWidth', 2, 'MarkerSize', 6);

grid off;
%title('Resolvent Singular Values vs Frequency \omega', 'FontSize', 22);
xlabel('Temporal frequency \omega', 'FontSize', 20);
ylabel('Singular value \sigma', 'FontSize', 20);

% legend({ ...
%     sprintf('kz = %.2f', kz_list(5)), ...
%     sprintf('kz = %.2f', kz_list(10)), ...
%     sprintf('kz = %.2f', kz_list(20)), ...
%     sprintf('kz = %.2f', kz_list(28)), ...
%     sprintf('kz = %.2f', kz_list(32))}, ...
%     'FontSize', 18, 'Location', 'northwest');

legend({ ...
    sprintf('kz = %.2f', kz_list(5)), ...
    sprintf('kz = %.2f', kz_list(10)), ...
    sprintf('kz = %.2f', kz_list(20)), ...
    sprintf('kz = %.2f', kz_list(28)), ...
    sprintf('kz = %.2f', kz_list(32))}, ...
    'FontSize', 18, ...
    'Orientation','horizontal', ...
    'NumColumns', 2, ...   % <-- first row: 3 items
    'Location','northwest');  % <-- inside, top-left



% Set axis tick label font size
set(gca, 'FontSize', 20);
lgd = legend;
set(lgd, 'Color', 'none');   % remove background
set(lgd, 'Box', 'off');      % optional: remove border
xlim([-1 1])
ylim([0 200])

%saveas(gcf, 'resolvent_spectrum_vs_omega.png');
filename = ['resolvent_spectrum_vs_omega_' num2str(Re) '.png'];
saveas(gcf, filename);

% --- Energy gain (sigma^2) vs phase speed c ---
figure;
plot(omega, sigma2_list_5 , 'o-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(omega, sigma2_list_10, 's--','LineWidth', 2, 'MarkerSize', 6);
plot(omega, sigma2_list_20, 'd-.','LineWidth', 2, 'MarkerSize', 6);
plot(omega, sigma2_list_28, 'd-.','LineWidth', 2, 'MarkerSize', 6);
plot(omega, sigma2_list_32, 'd-.', 'LineWidth', 2, 'MarkerSize', 6);

grid off;
%title('Resolvent Energy Gain vs Phase Speed', 'FontSize', 22);
xlabel('Temporal frequency \omega', 'FontSize', 20);
ylabel('Energy gain \sigma^{2}', 'FontSize', 20);

legend({ ...
    sprintf('kz = %.2f', kz_list(5)), ...
    sprintf('kz = %.2f', kz_list(10)), ...
    sprintf('kz = %.2f', kz_list(20)), ...
    sprintf('kz = %.2f', kz_list(28)), ...
    sprintf('kz = %.2f', kz_list(32))}, ...
    'FontSize', 18, 'Location', 'northwest');

% Set axis tick label font size
set(gca, 'FontSize', 20);
lgd = legend;
set(lgd, 'Color', 'none');   % remove background
set(lgd, 'Box', 'off');      % optional: remove border

filename = ['resolvent_energy_gain_kz_all_' num2str(Re) '.png'];
saveas(gcf, filename);

