%% === Confirm wall lies at mask midline (y where mask = 0.5) ===
% Make mask columns be along y:
M  = mask_smooth';               % Ny x Nx
Nx = numel(data_x); Ny = numel(data_y);

Mmax = max(M,[],1);               % per-column max (≈ mask_const)
Mn   = M ./ Mmax;                 % normalize to [0,1]

y50_bot = nan(1,Nx);
y50_top = nan(1,Nx);

for ix = 1:Nx
    mcol = Mn(:,ix);
    sgn  = mcol - 0.5;
    idxs = find(sgn(1:end-1).*sgn(2:end) <= 0);   % crossings of 0.5
    if numel(idxs) >= 2
        yhits = zeros(1,numel(idxs));
        for k = 1:numel(idxs)
            i  = idxs(k);
            y1 = data_y(i);   y2 = data_y(i+1);
            m1 = mcol(i);     m2 = mcol(i+1);
            yhits(k) = y1 + (0.5 - m1) * (y2 - y1) / (m2 - m1); % linear interp
        end
        yhits = sort(yhits);       % bottom, top
        y50_bot(ix) = yhits(1);
        y50_top(ix) = yhits(end);
    end
end

% Displacement length relative to analytic walls
ell_bot = y50_bot - y_wall2a(:).';
%ell_top = y50_top - y_wall2b(:).';

% Summaries (+ compare to your δ* = mask_threshold)
rms_bot = sqrt(nanmean(ell_bot.^2));  max_bot = nanmax(abs(ell_bot));
%rms_top = sqrt(nanmean(ell_top.^2));  max_top = nanmax(abs(ell_top));

fprintf('ℓ* bottom: RMS = %.3e (RMS/δ* = %.2f), max = %.3e (max/δ* = %.2f)\n', ...
    rms_bot, rms_bot/mask_threshold, max_bot, max_bot/mask_threshold);
% fprintf('ℓ* top   : RMS = %.3e (RMS/δ* = %.2f), max = %.3e (max/δ* = %.2f)\n', ...
%     rms_top, rms_top/mask_threshold, max_top, max_top/mask_threshold);

% Optional overlay plot
figure;
plot(data_x, y_wall2a, 'r-'); hold on;
plot(data_x, y50_bot,  'k.');
legend('wall^-','mask midline^-');
xlabel('x'); ylabel('y'); title('Wall vs mask midline (y_{50})');

% figure; hold on
% plot(data_x, ell_bot, 'b-', 'DisplayName','\ell^*_bottom');
% plot(data_x, ell_top, 'r-', 'DisplayName','\ell^*_top');
% yline( 0.1*mask_threshold,'k--','0.1\delta^*');
% yline(-0.1*mask_threshold,'k--','-0.1\delta^*');
% xlabel('x'); ylabel('\ell^*'); legend; grid on
% title('\ell^*(x) vs x');