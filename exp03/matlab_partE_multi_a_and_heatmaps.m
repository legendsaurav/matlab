%% ========================================================================
%  PART E -- MULTI-a WATERFALL/OVERSHOOT COMPARISON & MULTI-METRIC HEATMAPS
%  Reproduces report Figures 16-20 and Table 12 (report Sections 10.3-10.4).
%  This is NEW analysis added after the original Parts A-D delivery, in
%  response to extending the single-a=5 waterfall / single-metric heatmap
%  into multi-a and multi-metric comparisons.
%
%  Requires: Control System Toolbox (tf, step)
%  Needs on path: charEqCoeffs.m, dominantComplexPole.m, robustStepResponse.m,
%                 viridisLike.m
% ========================================================================
clc; clear; close all;

outdir = 'figures_partE';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ---- Figure 16: multi-a step-response waterfall (4 panels) -----------
a_compare = [2 5 9 20];
K_wf = logspace(log10(0.15), log10(300), 45);
Tend = 22.0; N = 2200;
t_common = linspace(0, Tend, N);

figure('visible', 'off', 'Position', [50 50 1000 850]);
cmapWF = viridisLike(numel(K_wf));
for p = 1:numel(a_compare)
    a = a_compare(p);
    subplot(2, 2, p);
    hold on;
    for i = 1:numel(K_wf)
        K = K_wf(i);
        sys = tf([K K], charEqCoeffs(a, K));
        y = step(sys, t_common);
        plot3(t_common, log10(K)*ones(size(t_common)), y, 'Color', cmapWF(i,:), 'LineWidth', 0.7);
    end
    xlabel('Time (s)', 'FontSize', 8); ylabel('log10(K)', 'FontSize', 8); zlabel('y(t)', 'FontSize', 8);
    title(sprintf('a = %d', a), 'FontWeight', 'bold');
    view(-60, 22); grid on; box on;
end
sgtitle('Continuous Step-Response Waterfall Across a  (K = 0.15 \rightarrow 300 in every panel)');
print(fullfile(outdir, 'fig16_waterfall_multi_a.png'), '-dpng', '-r150');

%% ---- Figure 17 + Table 12: true simulated PO(K) for several a --------
K_line = logspace(log10(0.15), log10(300), 70);
figure('visible', 'off', 'Position', [50 50 720 540]);
hold on; grid on; box on;
colorsA = [0.10 0.15 0.45; 0.30 0.30 0.40; 0.55 0.50 0.35; 0.85 0.72 0.20];
PO_at_top = zeros(size(a_compare));
for i = 1:numel(a_compare)
    a = a_compare(i);
    POs = zeros(size(K_line));
    for j = 1:numel(K_line)
        [~, y] = robustStepResponse(a, K_line(j));
        POs(j) = max((max(y) - 1) * 100, 0);
    end
    plot(K_line, POs, 'Color', colorsA(i,:), 'LineWidth', 2.1);
    PO_at_top(i) = POs(end);
end
set(gca, 'XScale', 'log');
xlabel('Gain, K (log scale)'); ylabel('True Simulated Peak Overshoot (%)');
title({'Simulated Overshoot vs. K for Several Values of a', '(same continuous K sweep as the waterfall panels)'});
legend(arrayfun(@(a) sprintf('a = %d', a), a_compare, 'UniformOutput', false), 'Location', 'northeast', 'FontSize', 9.5);
print(fullfile(outdir, 'fig17_po_vs_K_multi_a.png'), '-dpng', '-r150');

fprintf('\n=== Table 12: true simulated PO at the top of the sweep (K = %.1f) ===\n', K_line(end));
fprintf('%6s %20s\n', 'a', 'PO at K~300 (%)');
for i = 1:numel(a_compare)
    fprintf('%6d %20.1f\n', a_compare(i), PO_at_top(i));
end

%% ---- Figures 18-20: multi-metric heatmaps zeta / PO / Ts -------------
% One pass over the (K,a) grid computes the closed-loop pole ONCE per
% point; zeta, PO and Ts are all derived from that single pole, so all
% three heatmaps are computed from a single loop (not three separate
% sweeps) and are pixel-for-pixel consistent with each other.
K_grid = logspace(log10(0.4), log10(400), 140);
a_grid2 = linspace(1.0, 40, 140);

ZETA = zeros(numel(a_grid2), numel(K_grid));
PO   = zeros(size(ZETA));
TS   = zeros(size(ZETA));
for i = 1:numel(a_grid2)
    for j = 1:numel(K_grid)
        s = dominantComplexPole(a_grid2(i), K_grid(j));
        wn = abs(s); z = -real(s)/wn;
        ZETA(i,j) = z;
        if z < 1
            PO(i,j) = 100*exp(-z*pi/sqrt(1-z^2));
        else
            PO(i,j) = 0;
        end
        TS(i,j) = 3/max(z*wn, 1e-9);
    end
end

% --- Figure 18: zeta(K,a) ---
figure('visible', 'off', 'Position', [50 50 700 560]);
imagesc(log10(K_grid), a_grid2, ZETA);
set(gca, 'YDir', 'normal'); colormap(gca, 'parula');
cb = colorbar; ylabel(cb, 'Damping Ratio \zeta');
xlabel('log10(K)'); ylabel('Open-loop pole parameter, a');
title('Damping Ratio Sensitivity, \zeta(K, a)');
print(fullfile(outdir, 'fig18_heatmap_zeta.png'), '-dpng', '-r150');

% --- Figure 19: PO(K,a) ---
figure('visible', 'off', 'Position', [50 50 700 560]);
imagesc(log10(K_grid), a_grid2, PO);
set(gca, 'YDir', 'normal'); colormap(gca, 'hot');
cb = colorbar; ylabel(cb, 'Percent Overshoot (%)');
xlabel('log10(K)'); ylabel('Open-loop pole parameter, a');
title('Overshoot Sensitivity, PO(K, a)');
print(fullfile(outdir, 'fig19_heatmap_overshoot.png'), '-dpng', '-r150');

% --- Figure 20: Ts(K,a), clipped for display (true values kept in TS) ---
TS_display = min(TS, 20);   % Ts -> Inf as a->1 (zeta->0); clip for display only
figure('visible', 'off', 'Position', [50 50 700 560]);
imagesc(log10(K_grid), a_grid2, log10(TS_display));
set(gca, 'YDir', 'normal'); colormap(gca, 'bone');
cb = colorbar; ylabel(cb, 'log10( Settling Time / s ), clipped at 20 s');
xlabel('log10(K)'); ylabel('Open-loop pole parameter, a');
title('5% Settling-Time Sensitivity, Ts(K, a)');
print(fullfile(outdir, 'fig20_heatmap_settling_time.png'), '-dpng', '-r150');

%% ---- Verify the high-K asymptote Ts -> 6/(a-1) ------------------------
% Derived from the Sec. 6.2 high-gain pole asymptotics s ~ (1-a)/2 +- j*sqrt(K):
% the real part saturates at (1-a)/2 independent of K, so
% Ts = 3/(-Re(s)) -> 3/((a-1)/2) = 6/(a-1) as K -> Inf.
fprintf('\n=== Ts asymptote check: Ts -> 6/(a-1) as K -> Inf ===\n');
fprintf('%6s %14s %14s %14s\n', 'a', 'K', 'Ts (actual)', 'Ts (predicted)');
for a = [5 10 20]
    predicted = 6/(a-1);
    for K = [1000 10000 100000]
        s = dominantComplexPole(a, K);
        Ts_actual = 3/(-real(s));
        fprintf('%6d %14d %14.4f %14.4f\n', a, K, Ts_actual, predicted);
    end
end

fprintf('\nPart E complete. Figures saved to %s/\n', outdir);
