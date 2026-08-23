%% ========================================================================
%  PART C -- OVERSHOOT, SETTLING TIME, AND STEP-RESPONSE STUDIES
%  Reproduces report Figures 10, 11, 12, 13 and Tables 6, 7.
%
%  T(s) = K(s+1) / (s^3 + a*s^2 + K*s + K)
%
%  Requires: Control System Toolbox (tf, step)
%  Needs on path: charEqCoeffs.m, dominantComplexPole.m, realClosedLoopPole.m,
%                 viridisLike.m
% ========================================================================
clc; clear; close all;

outdir = 'figures_partC';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ---- Formula-based PO / Ts grid: a = 2..50, five gains ---------------
a_grid = 2:50;
K_values = [10000 20000 30000 50000 70000];

PO_formula = zeros(numel(a_grid), numel(K_values));
Ts_formula = zeros(size(PO_formula));
for j = 1:numel(K_values)
    K = K_values(j);
    for i = 1:numel(a_grid)
        a = a_grid(i);
        s = dominantComplexPole(a, K);
        wn = abs(s); zeta = -real(s)/wn;
        PO_formula(i,j) = 100*exp(-zeta*pi/sqrt(1-zeta^2));
        Ts_formula(i,j) = 3/(zeta*wn);
    end
end

%% ---- Figure 10: PO and Ts vs a (formula), five gains -----------------
cmap = viridisLike(numel(K_values));
figure('visible', 'off', 'Position', [50 50 1250 480]);

subplot(1,2,1); hold on; grid on; box on;
for j = 1:numel(K_values)
    plot(a_grid, PO_formula(:,j), 'Color', cmap(j,:), 'LineWidth', 2);
end
xlabel('Open-Loop Pole Position, a'); ylabel('Peak Overshoot (%)');
title('Peak Overshoot vs. a (2nd-order formula)');
legend(arrayfun(@(K) sprintf('K = %d', K), K_values, 'UniformOutput', false), 'Location', 'northeast', 'FontSize', 9);

subplot(1,2,2); hold on; grid on; box on;
for j = 1:numel(K_values)
    plot(a_grid, Ts_formula(:,j), 'Color', cmap(j,:), 'LineWidth', 2);
end
xlabel('Open-Loop Pole Position, a'); ylabel('5%% Settling Time (s)');
title('5%% Settling Time vs. a (2nd-order formula)');
legend(arrayfun(@(K) sprintf('K = %d', K), K_values, 'UniformOutput', false), 'Location', 'northeast', 'FontSize', 9);
print(fullfile(outdir, 'fig10_overshoot_settling_vs_a.png'), '-dpng', '-r150');

%% ---- Benchmark against true simulated step response ------------------
% This is the "core" computation from the report's Section 8.2 code
% block, looped over the full (a,K) grid to build Table 6/7 and Figure 11.
% Uses robustStepResponse() (explicit fine time grid) rather than handing
% step() only a scalar final time -- see that function's header comment
% for why: an under-resolved grid can miss a sharp, lightly-damped first
% peak by close to a percentage point.
PO_actual = zeros(size(PO_formula));
Ts_actual = zeros(size(PO_formula));
s3_mean_track = zeros(numel(K_values), 1);

fprintf('\n=== Table 6 (subset): formula vs. simulated overshoot/settling ===\n');
fprintf('%4s %8s %10s %10s %10s %10s\n', 'a', 'K', 'PO_form', 'PO_sim', 'Ts_form', 'Ts_sim');
for j = 1:numel(K_values)
    K = K_values(j);
    s3_vals = zeros(numel(a_grid),1);
    for i = 1:numel(a_grid)
        a = a_grid(i);

        % ---- True (non-approximated) step response ----
        [t, y] = robustStepResponse(a, K);

        PO_actual(i,j) = max((max(y)-1)*100, 0);
        band = 0.05;
        outside = find(abs(y-1) > band);
        if isempty(outside), Ts_actual(i,j) = 0; else, Ts_actual(i,j) = t(outside(end)); end

        s3_vals(i) = realClosedLoopPole(a, K);

        if ismember(a, [2 9 20 50])
            fprintf('%4d %8d %10.2f %10.2f %10.4f %10.4f\n', a, K, ...
                PO_formula(i,j), PO_actual(i,j), Ts_formula(i,j), Ts_actual(i,j));
        end
    end
    s3_mean_track(j) = mean(s3_vals);
end

fprintf('\n=== Table 7: mean approximation error vs. K ===\n');
fprintf('%8s %14s %14s %12s\n', 'K', 'mean|POerr|', 'mean|Tserr|', 'mean s3');
for j = 1:numel(K_values)
    po_err = mean(abs(PO_formula(:,j) - PO_actual(:,j)));
    ts_err = mean(abs(Ts_formula(:,j) - Ts_actual(:,j)));
    fprintf('%8d %14.3f %14.4f %12.4f\n', K_values(j), po_err, ts_err, s3_mean_track(j));
end

%% ---- Figure 11: parity plots, formula vs. true simulated -------------
figure('visible', 'off', 'Position', [50 50 1250 520]);
subplot(1,2,1); hold on; grid on; box on;
for j = 1:numel(K_values)
    scatter(PO_formula(:,j), PO_actual(:,j), 18, cmap(j,:), 'filled');
end
lims = [0 100];
plot(lims, lims, 'k--', 'LineWidth', 1);
xlim(lims); ylim(lims);
xlabel('Peak Overshoot -- formula (%)'); ylabel('Peak Overshoot -- simulated (%)');
title('Overshoot: Formula vs. True Simulated Response');

subplot(1,2,2); hold on; grid on; box on;
for j = 1:numel(K_values)
    scatter(Ts_formula(:,j), Ts_actual(:,j), 18, cmap(j,:), 'filled');
end
lims2 = [0 max(Ts_formula(:))*1.05];
plot(lims2, lims2, 'k--', 'LineWidth', 1);
xlim(lims2); ylim(lims2);
xlabel('Settling Time -- formula (s)'); ylabel('Settling Time -- simulated (s)');
title('Settling Time: Formula vs. True Simulated Response');
print(fullfile(outdir, 'fig11_formula_vs_actual_parity.png'), '-dpng', '-r150');

%% ---- Figure 12: third closed-loop pole s3 -> -1 as K -> Inf ----------
a_fix_list = [2 5 9 20];
K_sweep = logspace(1, log10(2e5), 60);
figure('visible', 'off', 'Position', [50 50 720 520]);
hold on; grid on; box on;
cmap3 = [0.10 0.15 0.45; 0.35 0.35 0.45; 0.55 0.50 0.35; 0.85 0.72 0.20];
for i = 1:numel(a_fix_list)
    a = a_fix_list(i);
    s3 = arrayfun(@(K) realClosedLoopPole(a, K), K_sweep);
    plot(K_sweep, s3, 'Color', cmap3(i,:), 'LineWidth', 2);
end
yline(-1, 'r--', 'LineWidth', 1.3);
set(gca, 'XScale', 'log');
xlabel('Gain, K (log scale)'); ylabel('Real closed-loop pole  s_3');
title({'Third Closed-Loop Pole Approaches the Zero as K \rightarrow \infty', ...
       '(pole/zero near-cancellation justifies the 2nd-order approximation)'});
legend([arrayfun(@(a) sprintf('a = %d', a), a_fix_list, 'UniformOutput', false), {'open-loop zero, s=-1'}], ...
    'Location', 'southeast', 'FontSize', 9);
print(fullfile(outdir, 'fig12_s3_approaches_zero.png'), '-dpng', '-r150');

%% ---- Figure 13: step response, low gain (K=5) vs high gain (K=10000) -
a_case = [2 5 9 10];
colors2 = [0.23 0.44 0.63; 0.76 0.27 0.05; 0.25 0.62 0.30; 0.54 0.31 0.75];

figure('visible', 'off', 'Position', [50 50 1300 500]);
subplot(1,2,1); hold on; grid on; box on;
for i = 1:numel(a_case)
    [t, y] = robustStepResponse(a_case(i), 5, 14);
    plot(t, y, 'Color', colors2(i,:), 'LineWidth', 1.6);
end
yline(1, 'k:', 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('Output, y(t)');
title('Step Response -- Low Gain, K = 5');
legend(arrayfun(@(a) sprintf('a = %d', a), a_case, 'UniformOutput', false), 'Location', 'best', 'FontSize', 9);

subplot(1,2,2); hold on; grid on; box on;
for i = 1:numel(a_case)
    [t, y] = robustStepResponse(a_case(i), 10000, 3);
    plot(t, y, 'Color', colors2(i,:), 'LineWidth', 1.1);
end
yline(1, 'k:', 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('Output, y(t)');
title('Step Response -- High Gain, K = 10,000');
legend(arrayfun(@(a) sprintf('a = %d', a), a_case, 'UniformOutput', false), 'Location', 'best', 'FontSize', 9);
print(fullfile(outdir, 'fig13_step_response_low_high_gain.png'), '-dpng', '-r150');

fprintf('\nPart C complete. Figures saved to %s/\n', outdir);
