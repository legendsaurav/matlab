%% ========================================================================
%  PART D -- EXTENDED ANALYSIS: OPTIMAL DAMPING GAIN & RINGING
%  Reproduces report Figures 14, 15 and Table 8.
%
%  NOTE: this script previously also generated a single-a=5 waterfall and
%  a single-metric PO(K,a) heatmap (old Figures 16-17). Those have been
%  superseded by the multi-a and multi-metric versions in
%  matlab_partE_multi_a_and_heatmaps.m and were removed from here to
%  avoid two different, conflicting "Figure 16"/"Figure 17" outputs.
%
%  Requires: Control System Toolbox (tf, step)
%  Needs on path: charEqCoeffs.m, dominantComplexPole.m, magmaLike.m,
%                 robustStepResponse.m
% ========================================================================
clc; clear; close all;

outdir = 'figures_partD';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ---- Figure 14: zeta(K) is non-monotonic; optimal gain K* -----------
% Restricted to a < 9 (see report Sec. 10.1): for a > 9 there is a
% window of K with three REAL closed-loop poles (the loop of Part A/
% Section 5.3), where the complex-pole damping ratio isn't meaningful.
a_list = [2 3 5 7];
cmapD = magmaLike(numel(a_list));

figure('visible', 'off', 'Position', [50 50 760 540]);
hold on; grid on; box on;
Kstar = zeros(size(a_list));
zetaMax = zeros(size(a_list));
for i = 1:numel(a_list)
    a = a_list(i);
    Ks = logspace(-3, 5, 400);
    zetas = arrayfun(@(K) -real(dominantComplexPole(a,K))/abs(dominantComplexPole(a,K)), Ks);
    plot(Ks, zetas, 'Color', cmapD(i,:), 'LineWidth', 2);

    % bounded 1-D search for the interior maximum, in log(K)
    negZeta = @(lk) -(-real(dominantComplexPole(a, exp(lk))) / abs(dominantComplexPole(a, exp(lk))));
    [lk_opt, negz] = fminbnd(negZeta, log(1e-3), log(1e5));
    Kstar(i) = exp(lk_opt);
    zetaMax(i) = -negz;
    plot(Kstar(i), zetaMax(i), 'o', 'Color', cmapD(i,:), 'MarkerFaceColor', cmapD(i,:), ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 7);
end
set(gca, 'XScale', 'log');
xlabel('Gain, K (log scale)'); ylabel('Damping Ratio, \zeta');
title('\zeta vs. K is Non-Monotonic -- an Optimal-Damping Gain K^* Exists');
legend(arrayfun(@(a) sprintf('a = %d', a), a_list, 'UniformOutput', false), 'Location', 'northeast', 'FontSize', 9.5);
print(fullfile(outdir, 'fig14_zeta_vs_K_nonmonotonic.png'), '-dpng', '-r150');

fprintf('\n=== Table 8: best-damping gain K* vs. heuristic (a-1)*sqrt(a) ===\n');
fprintf('%6s %14s %12s %18s\n', 'a', 'K*', 'zeta_max', 'heuristic(a-1)sqrt(a)');
for i = 1:numel(a_list)
    a = a_list(i);
    fprintf('%6d %14.3f %12.4f %18.3f\n', a, Kstar(i), zetaMax(i), (a-1)*sqrt(a));
end

%% ---- Figure 15: overshoot with/without ringing, a=12 -----------------
a_demo = 12;
K_demo = [41 200];
demo_labels = {'K = 41  (inside loop: 3 real poles \rightarrow single-hump overshoot)', ...
               'K = 200 (past loop: complex poles \rightarrow overshoot WITH ringing)'};
figure('visible', 'off', 'Position', [50 50 760 540]);
hold on; grid on; box on;
colorsR = [0.76 0.27 0.05; 0.25 0.62 0.30];
hh = gobjects(1,2);
fprintf('\n=== Ringing demo (a=12): poles, overshoot, ring-crossings ===\n');
for i = 1:2
    K = K_demo(i);
    [t, y] = robustStepResponse(a_demo, K, 3);
    hh(i) = plot(t, y, 'Color', colorsR(i,:), 'LineWidth', 1.8);
    r = roots(charEqCoeffs(a_demo, K));
    all_real = all(abs(imag(r)) < 1e-6);
    [~, peak_idx] = max(y);
    rings = sum(diff(sign(y(peak_idx:end) - 1)) ~= 0);
    fprintf('K=%3d  poles=[%s]  all_real=%d  PO=%.2f%%  ring_crossings=%d\n', ...
        K, sprintf('%.3f%+.3fi  ', [real(r) imag(r)]'), all_real, (max(y)-1)*100, rings);
end
yline(1, 'k:', 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('Output, y(t)');
title(sprintf('Overshoot Without Ringing vs. With Ringing  (a = %d)', a_demo));
legend(hh, demo_labels, 'Location', 'best', 'FontSize', 8.5);
print(fullfile(outdir, 'fig15_ringing_vs_nonringing.png'), '-dpng', '-r150');


fprintf('\nPart D complete. Figures saved to %s/\n', outdir);
