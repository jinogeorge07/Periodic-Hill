% Amplification factor and growth
clc; clear all; close all;

folderpath = './'

Re = [50,100,190];

c_number = 12;
compute_deviation = "false"
%% Wavenumbers/frequency
kxn = 1; kzn = 36; kx = 1;
omega = 1;
kx_list = logspace(-4,0.48,kxn);
kz_list = logspace(-2,1.2,kzn);
c_list  = linspace(-1,1,c_number);
omega = -c_list * kx;

for i = 1:length(Re)
 
filename = ['sigma_list_Re' num2str(Re(i)) '_kz5.mat'];
data5 = load(filename);
filename = ['sigma_list_Re' num2str(Re(i)) '_kz10.mat'];
data10 = load(filename);
filename = ['sigma_list_Re' num2str(Re(i)) '_kz20.mat'];
data20 = load(filename);
filename = ['sigma_list_Re' num2str(Re(i)) '_kz28.mat'];
data28 = load(filename);
filename = ['sigma_list_Re' num2str(Re(i)) '_kz32.mat'];
data32 = load(filename);

writematrix(data5.sigma_list,['sigma_list_Re' num2str(Re(i)) '_kz5.csv'])
writematrix(data10.sigma_list,['sigma_list_Re' num2str(Re(i)) '_kz10.csv'])
writematrix(data20.sigma_list,['sigma_list_Re' num2str(Re(i)) '_kz20.csv'])
writematrix(data28.sigma_list,['sigma_list_Re' num2str(Re(i)) '_kz28.csv'])
writematrix(data32.sigma_list,['sigma_list_Re' num2str(Re(i)) '_kz32.csv'])

sigma_list_5 = data5.sigma_list;
sigma_list_10 = data10.sigma_list;
sigma_list_20 = data20.sigma_list;
sigma_list_28 = data28.sigma_list;
sigma_list_32 = data32.sigma_list;

sigma2_list_5  = sigma_list_5.^2;
sigma2_list_10 = sigma_list_10.^2;
sigma2_list_20 = sigma_list_20.^2;
sigma2_list_28 = sigma_list_28.^2;
sigma2_list_32 = sigma_list_32.^2;

% Amplification gain------------------------------------------------------

% Set axis tick label font size
set(gca, 'FontSize', 20);
lgd = legend;
set(lgd, 'Color', 'none');   % remove background
set(lgd, 'Box', 'off');      % optional: remove border
xlim([-1 1])
ylim([0 200])

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
filename = ['resolvent_spectrum_vs_omega_' num2str(Re(i)) '.png'];
saveas(gcf, filename);

end
