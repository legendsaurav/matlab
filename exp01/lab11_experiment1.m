%% EE208 Control Systems Lab Experiment 1
% Open-loop transfer function:
% G(s) = K(1+s)(1+5s) / [s^2(s^2 + 2s + 10)]
%
% This script performs root locus, continuous-gain, step-response,
% second-order approximation, transition, and pole-variation analyses.
% All computed data is exported with writetable().

clear; close all; clc;

outputDir = pwd;
baseKs = [0.1, 0.5, 1, 2, 5, 10, 20, 50];
continuousKs = unique([logspace(log10(0.01), log10(300), 240), baseKs, 2.36]);
transitionKs = linspace(1, 5, 200);
baseOL = tf([5 6 1], [1 2 10 0 0]);

% -------------------------------------------------------------------------
% ROOT LOCUS ANALYSIS
% -------------------------------------------------------------------------
rootRows = cell(numel(baseKs) * 4, 4);
rowIdx = 1;
for kIdx = 1:numel(baseKs)
    K = baseKs(kIdx);
    [polesSorted, ~, ~, ~, ~, ~, ~] = analyzeClosedLoop(K, 'nominal');
    for pIdx = 1:4
        rootRows{rowIdx, 1} = K;
        rootRows{rowIdx, 2} = pIdx;
        rootRows{rowIdx, 3} = real(polesSorted(pIdx));
        rootRows{rowIdx, 4} = imag(polesSorted(pIdx));
        rowIdx = rowIdx + 1;
    end
end
rootTable = cell2table(rootRows, 'VariableNames', {'K', 'PoleIndex', 'Pole_Real', 'Pole_Imag'});
writetable(rootTable, fullfile(outputDir, 'root_locus_poles.csv'));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1450 650]);
subplot(1, 2, 1);
rlocus(baseOL);
hold on;
plot(real(pole(baseOL)), imag(pole(baseOL)), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(zero(baseOL)), imag(zero(baseOL)), 'go', 'MarkerSize', 10, 'LineWidth', 2);
for kIdx = 1:numel(baseKs)
    K = baseKs(kIdx);
    polesSorted = analyzeClosedLoop(K, 'nominal');
    plot(real(polesSorted), imag(polesSorted), 'k.', 'MarkerSize', 18);
    text(real(polesSorted(1)) + 0.05, imag(polesSorted(1)), sprintf('  K=%.1f', K), 'FontSize', 8);
end
grid on;
title('Root Locus of G(s) with Sample Closed-Loop Poles');
xlabel('Real Axis'); ylabel('Imaginary Axis');
legend('Root locus', 'Open-loop poles', 'Open-loop zeros', 'Sample closed-loop poles', 'Location', 'best');

subplot(1, 2, 2);
rlocus(baseOL);
hold on;
plot(real(pole(baseOL)), imag(pole(baseOL)), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(zero(baseOL)), imag(zero(baseOL)), 'go', 'MarkerSize', 10, 'LineWidth', 2);
for kIdx = 1:numel(baseKs)
    K = baseKs(kIdx);
    polesSorted = analyzeClosedLoop(K, 'nominal');
    plot(real(polesSorted), imag(polesSorted), 'k.', 'MarkerSize', 18);
end
grid on;
xlim([-2.5 0.5]);
ylim([-6 6]);
title('Zoomed Root Locus Near the Origin');
xlabel('Real Axis'); ylabel('Imaginary Axis');
saveas(fig, fullfile(outputDir, 'fig_root_locus.png'));
close(fig);

% -------------------------------------------------------------------------
% CONTINUOUS GAIN ANALYSIS
% -------------------------------------------------------------------------
continuousRows = analyzeGainSweep(continuousKs, 'nominal');
continuousTable = struct2table(continuousRows);
writetable(continuousTable, fullfile(outputDir, 'gain_analysis_continuous.csv'));

% Generate the pole-variation and comprehensive outputs early so they are
% not dependent on the long tail of this script.
run(fullfile(outputDir, 'lab11_experiment1_tail.m'));

transitionRows = analyzeGainSweep(transitionKs, 'nominal');
transitionTable = struct2table(transitionRows);
writetable(transitionTable, fullfile(outputDir, 'transition_analysis.csv'));

% -------------------------------------------------------------------------
% STEP RESPONSE SIMULATIONS FOR SPECIFIC K VALUES
% -------------------------------------------------------------------------
measuredRows = cell(numel(baseKs), 6);
approxRows = cell(numel(baseKs), 11);
comparisonRows = cell(numel(baseKs), 14);
stepSeries = struct();

for kIdx = 1:numel(baseKs)
    K = baseKs(kIdx);
    [polesSorted, dominantPole, zeta, omegaN, omegaD, numComplexPairs, ~] = analyzeClosedLoop(K, 'nominal');
    [t, y, metrics] = simulateAndMeasureStep(K, 'nominal', polesSorted, dominantPole);

    measuredRows{kIdx, 1} = K;
    measuredRows{kIdx, 2} = metrics.PO;
    measuredRows{kIdx, 3} = metrics.tp;
    measuredRows{kIdx, 4} = metrics.ts;
    measuredRows{kIdx, 5} = metrics.ess;
    measuredRows{kIdx, 6} = metrics.yss;

    approxRows{kIdx, 1} = K;
    approxRows{kIdx, 2} = zeta;
    approxRows{kIdx, 3} = omegaN;
    approxRows{kIdx, 4} = omegaD;
    approxRows{kIdx, 5} = metrics.PO_approx;
    approxRows{kIdx, 6} = metrics.tp_approx;
    approxRows{kIdx, 7} = metrics.ts_approx;
    approxRows{kIdx, 8} = metrics.ess_approx;
    approxRows{kIdx, 9} = numComplexPairs;
    approxRows{kIdx, 10} = real(dominantPole);
    approxRows{kIdx, 11} = imag(dominantPole);

    comparisonRows{kIdx, 1} = K;
    comparisonRows{kIdx, 2} = zeta;
    comparisonRows{kIdx, 3} = omegaN;
    comparisonRows{kIdx, 4} = metrics.PO;
    comparisonRows{kIdx, 5} = metrics.PO_approx;
    comparisonRows{kIdx, 6} = safePercentError(metrics.PO, metrics.PO_approx);
    comparisonRows{kIdx, 7} = metrics.tp;
    comparisonRows{kIdx, 8} = metrics.tp_approx;
    comparisonRows{kIdx, 9} = safePercentError(metrics.tp, metrics.tp_approx);
    comparisonRows{kIdx, 10} = metrics.ts;
    comparisonRows{kIdx, 11} = metrics.ts_approx;
    comparisonRows{kIdx, 12} = safePercentError(metrics.ts, metrics.ts_approx);
    comparisonRows{kIdx, 13} = metrics.ess;
    comparisonRows{kIdx, 14} = metrics.ess_approx;

    stepSeries.(sprintf('K_%s', sanitizeValue(K))) = struct('t', t, 'y', y);
end

measuredTable = cell2table(measuredRows, 'VariableNames', {'K', 'PO_Measured', 'tp_Measured', 'ts_Measured', 'ess_Measured', 'y_ss'});
writetable(measuredTable, fullfile(outputDir, 'measured_metrics.csv'));

approxTable = cell2table(approxRows, 'VariableNames', {'K', 'zeta', 'omega_n', 'omega_d', 'PO_Approx', 'tp_Approx', 'ts_Approx', 'ess_Approx', 'NumComplexPairs', 'DominantPole_Real', 'DominantPole_Imag'});
writetable(approxTable, fullfile(outputDir, 'second_order_approx.csv'));

comparisonTable = cell2table(comparisonRows, 'VariableNames', {'K', 'zeta', 'omega_n', 'PO_Measured', 'PO_Approx', 'PO_Error_Percent', 'tp_Measured', 'tp_Approx', 'tp_Error_Percent', 'ts_Measured', 'ts_Approx', 'ts_Error_Percent', 'ess_Measured', 'ess_Approx'});
writetable(comparisonTable, fullfile(outputDir, 'comparison_table.csv'));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 900]);
for kIdx = 1:numel(baseKs)
    subplot(2, 4, kIdx);
    series = stepSeries.(sprintf('K_%s', sanitizeValue(baseKs(kIdx))));
    plot(series.t, series.y, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)'); ylabel('y(t)');
    title(sprintf('K = %.1f', baseKs(kIdx)));
end
sgtitle('Unit Step Responses for the Eight Selected Gains');
saveas(fig, fullfile(outputDir, 'fig_step_responses.png'));
close(fig);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 650]);
subplot(1, 3, 1);
semilogx(baseKs, cell2mat(measuredRows(:, 2)), 'bo-', 'LineWidth', 1.5); hold on;
semilogx(baseKs, cell2mat(comparisonRows(:, 5)), 'rs--', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('PO (%)'); title('Percentage Overshoot'); legend('Measured', '2nd-Order Approx', 'Location', 'best');

subplot(1, 3, 2);
semilogx(baseKs, cell2mat(measuredRows(:, 3)), 'bo-', 'LineWidth', 1.5); hold on;
semilogx(baseKs, cell2mat(comparisonRows(:, 8)), 'rs--', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('t_p (s)'); title('Peak Time'); legend('Measured', '2nd-Order Approx', 'Location', 'best');

subplot(1, 3, 3);
semilogx(baseKs, cell2mat(measuredRows(:, 4)), 'bo-', 'LineWidth', 1.5); hold on;
semilogx(baseKs, cell2mat(comparisonRows(:, 11)), 'rs--', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('t_s (s)'); title('Settling Time'); legend('Measured', '2nd-Order Approx', 'Location', 'best');
sgtitle('Measured Versus 2nd-Order Approximation');
saveas(fig, fullfile(outputDir, 'fig_comparison.png'));
close(fig);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 650]);
subplot(1, 3, 1);
semilogx(baseKs, cell2mat(comparisonRows(:, 6)), 'm-o', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('Error (%)'); title('PO Error');

subplot(1, 3, 2);
semilogx(baseKs, cell2mat(comparisonRows(:, 9)), 'm-o', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('Error (%)'); title('Peak Time Error');

subplot(1, 3, 3);
semilogx(baseKs, cell2mat(comparisonRows(:, 12)), 'm-o', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('Error (%)'); title('Settling Time Error');
sgtitle('Approximation Error Trends');
saveas(fig, fullfile(outputDir, 'fig_error_analysis.png'));
close(fig);

% -------------------------------------------------------------------------
% TRANSITION ANALYSIS PLOTS
% -------------------------------------------------------------------------
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1550 900]);

subplot(2, 3, 1);
plot(transitionTable.K, transitionTable.num_complex_pairs, 'k-', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('Pairs'); title('Number of Complex Pairs'); xline(2.36, 'r--', 'LineWidth', 1.2);

subplot(2, 3, 2);
plot(transitionTable.K, transitionTable.zeta, 'b-', 'LineWidth', 1.5); hold on;
plot(transitionTable.K, transitionTable.omega_n, 'r--', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('Value'); title('zeta and omega_n'); legend('zeta', 'omega_n', 'Location', 'best'); xline(2.36, 'k:', 'LineWidth', 1.2);

subplot(2, 3, 3);
plot(transitionTable.K, transitionTable.PO_measured, 'b-', 'LineWidth', 1.5); hold on;
plot(transitionTable.K, transitionTable.PO_approx, 'r--', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('PO (%)'); title('Overshoot Across Transition'); legend('Measured', 'Approx', 'Location', 'best'); xline(2.36, 'k:', 'LineWidth', 1.2);

subplot(2, 3, 4);
plot(transitionTable.K, transitionTable.tp_measured, 'b-', 'LineWidth', 1.5); hold on;
plot(transitionTable.K, transitionTable.tp_approx, 'r--', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('t_p (s)'); title('Peak Time Across Transition'); legend('Measured', 'Approx', 'Location', 'best'); xline(2.36, 'k:', 'LineWidth', 1.2);

subplot(2, 3, 5);
plot(transitionTable.K, transitionTable.ts_measured, 'b-', 'LineWidth', 1.5); hold on;
plot(transitionTable.K, transitionTable.ts_approx, 'r--', 'LineWidth', 1.5);
grid on; xlabel('K'); ylabel('t_s (s)'); title('Settling Time Across Transition'); legend('Measured', 'Approx', 'Location', 'best'); xline(2.36, 'k:', 'LineWidth', 1.2);

subplot(2, 3, 6);
plot(transitionTable.K, transitionTable.P1_real, 'LineWidth', 1.2); hold on;
plot(transitionTable.K, transitionTable.P2_real, 'LineWidth', 1.2);
plot(transitionTable.K, transitionTable.P3_real, 'LineWidth', 1.2);
plot(transitionTable.K, transitionTable.P4_real, 'LineWidth', 1.2);
plot(transitionTable.K, transitionTable.P1_imag, '--', 'LineWidth', 1.2);
plot(transitionTable.K, transitionTable.P2_imag, '--', 'LineWidth', 1.2);
plot(transitionTable.K, transitionTable.P3_imag, '--', 'LineWidth', 1.2);
plot(transitionTable.K, transitionTable.P4_imag, '--', 'LineWidth', 1.2);
grid on; xlabel('K'); ylabel('Pole Parts'); title('Pole Trajectories'); xline(2.36, 'k:', 'LineWidth', 1.2);

sgtitle('Focused Transition Analysis from K = 1 to 5');
saveas(fig, fullfile(outputDir, 'fig_transition_analysis.png'));
close(fig);

% -------------------------------------------------------------------------
% POLE VARIATION ANALYSIS AT K = 5
% -------------------------------------------------------------------------
poleVariations = {
    struct('label', 'Negative20', 'a', 1.6, 'b', 6.4),
    struct('label', 'Nominal',    'a', 2.0, 'b', 10.0),
    struct('label', 'Positive20', 'a', 2.4, 'b', 12.96)
};

variationRows = cell(3, 14);
variationSeries = struct();
variationPoles = struct();
variationDebugFile = fullfile(outputDir, 'pole_variation_debug.txt');
variationDebugId = fopen(variationDebugFile, 'w');

try
    fprintf(variationDebugId, 'Pole variation block start\n');
    for vIdx = 1:numel(poleVariations)
        variation = poleVariations{vIdx};
        fprintf(variationDebugId, 'Begin %s\n', variation.label);

        [polesSorted, dominantPole, ~, ~, ~, ~, ~] = analyzeClosedLoop(5, variation.label);
        [t, y, metrics] = simulateAndMeasureStep(5, variation.label, polesSorted, dominantPole);
        fprintf(variationDebugId, 'Simulated %s\n', variation.label);

        variationRows{vIdx, 1} = variation.label;
        variationRows{vIdx, 2} = 5;
        variationRows{vIdx, 3} = metrics.PO;
        variationRows{vIdx, 4} = metrics.tp;
        variationRows{vIdx, 5} = metrics.ts;
        variationRows{vIdx, 6} = metrics.ess;
        variationRows{vIdx, 7} = real(polesSorted(1));
        variationRows{vIdx, 8} = imag(polesSorted(1));
        variationRows{vIdx, 9} = real(polesSorted(2));
        variationRows{vIdx, 10} = imag(polesSorted(2));
        variationRows{vIdx, 11} = real(polesSorted(3));
        variationRows{vIdx, 12} = imag(polesSorted(3));
        variationRows{vIdx, 13} = real(polesSorted(4));
        variationRows{vIdx, 14} = imag(polesSorted(4));

        fieldPrefix = sanitizeText(variation.label);
        variationSeries.([fieldPrefix '_t']) = t;
        variationSeries.([fieldPrefix '_y']) = y;
        variationPoles.(variation.label) = polesSorted;
    end

    fprintf(variationDebugId, 'Loop complete\n');
    variationTable = cell2table(variationRows, 'VariableNames', {'Variation_Label', 'K', 'PO', 'tp', 'ts', 'ess', 'Pole1_Real', 'Pole1_Imag', 'Pole2_Real', 'Pole2_Imag', 'Pole3_Real', 'Pole3_Imag', 'Pole4_Real', 'Pole4_Imag'});
    writetable(variationTable, fullfile(outputDir, 'pole_variation.csv'));

    timeGrid = buildCommonTimeGrid({variationSeries.Negative20_t, variationSeries.Nominal_t, variationSeries.Positive20_t});
    respNeg20 = interp1(variationSeries.Negative20_t, variationSeries.Negative20_y, timeGrid, 'linear', 'extrap');
    respNominal = interp1(variationSeries.Nominal_t, variationSeries.Nominal_y, timeGrid, 'linear', 'extrap');
    respPos20 = interp1(variationSeries.Positive20_t, variationSeries.Positive20_y, timeGrid, 'linear', 'extrap');

    timeseriesTable = table(timeGrid(:), respNeg20(:), respNominal(:), respPos20(:), 'VariableNames', {'Time', 'Response_Neg20', 'Response_Nominal', 'Response_Pos20'});
    writetable(timeseriesTable, fullfile(outputDir, 'pole_variation_timeseries.csv'));

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 600]);
    plot(timeGrid, respNeg20, 'r-', 'LineWidth', 1.5); hold on;
    plot(timeGrid, respNominal, 'k--', 'LineWidth', 1.5);
    plot(timeGrid, respPos20, 'b-.', 'LineWidth', 1.5);
    grid on; xlabel('Time (s)'); ylabel('y(t)'); title('Pole Variation Step Response Comparison'); legend('Negative 20%', 'Nominal', 'Positive 20%', 'Location', 'best');
    saveas(fig, fullfile(outputDir, 'fig_pole_variation_step.png'));
    close(fig);

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 850]);
    colors = lines(3);
    for vIdx = 1:numel(poleVariations)
        variation = poleVariations{vIdx};
        polesSorted = variationPoles.(variation.label);
        plot(real(polesSorted), imag(polesSorted), 'o', 'Color', colors(vIdx, :), 'MarkerSize', 9, 'LineWidth', 1.8); hold on;
        for pIdx = 1:numel(polesSorted)
            text(real(polesSorted(pIdx)) + 0.05, imag(polesSorted(pIdx)), sprintf(' %s', variation.label), 'Color', colors(vIdx, :));
        end
    end
    grid on; xlabel('Real Axis'); ylabel('Imaginary Axis'); title('Closed-Loop Pole Locations for Pole Variations at K = 5');
    saveas(fig, fullfile(outputDir, 'fig_pole_variation_poles.png'));
    close(fig);

    fprintf(variationDebugId, 'Pole variation block complete\n');
catch ME
    if variationDebugId > 0
        fprintf(variationDebugId, 'ERROR: %s\n', ME.message);
        fprintf(variationDebugId, '%s\n', getReport(ME, 'extended'));
    end
    if variationDebugId > 0
        fclose(variationDebugId);
    end
    error('Pole variation analysis failed: %s', ME.message);
end

if variationDebugId > 0
    fclose(variationDebugId);
end

% -------------------------------------------------------------------------
% COMPREHENSIVE GAIN ANALYSIS FIGURE
% -------------------------------------------------------------------------
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 1700 1300]);

subplot(4, 3, 1);
plot(continuousTable.K, continuousTable.PO_measured, 'b-', 'LineWidth', 1.2); grid on; xlabel('K'); ylabel('PO (%)'); title('PO vs K'); xline(2.36, 'k:');

subplot(4, 3, 2);
plot(continuousTable.K, continuousTable.tp_measured, 'b-', 'LineWidth', 1.2); grid on; xlabel('K'); ylabel('t_p (s)'); title('tp vs K'); xline(2.36, 'k:');

subplot(4, 3, 3);
plot(continuousTable.K, continuousTable.ts_measured, 'b-', 'LineWidth', 1.2); grid on; xlabel('K'); ylabel('t_s (s)'); title('ts vs K'); xline(2.36, 'k:');

subplot(4, 3, 4);
plot(continuousTable.K, continuousTable.zeta, 'r-', 'LineWidth', 1.2); grid on; xlabel('K'); ylabel('zeta'); title('Damping Ratio'); xline(2.36, 'k:');

subplot(4, 3, 5);
plot(continuousTable.K, continuousTable.omega_n, 'r-', 'LineWidth', 1.2); grid on; xlabel('K'); ylabel('\omega_n'); title('Natural Frequency'); xline(2.36, 'k:');

subplot(4, 3, 6);
plot(continuousTable.K, continuousTable.omega_d, 'r-', 'LineWidth', 1.2); grid on; xlabel('K'); ylabel('\omega_d'); title('Damped Frequency'); xline(2.36, 'k:');

subplot(4, 3, 7);
plot(continuousTable.K, continuousTable.P1_real, 'LineWidth', 1.2); hold on;
plot(continuousTable.K, continuousTable.P2_real, 'LineWidth', 1.2);
plot(continuousTable.K, continuousTable.P3_real, 'LineWidth', 1.2);
plot(continuousTable.K, continuousTable.P4_real, 'LineWidth', 1.2);
grid on; xlabel('K'); ylabel('Real Part'); title('Pole Real Parts'); xline(2.36, 'k:');

subplot(4, 3, 8);
plot(continuousTable.K, continuousTable.P1_imag, 'LineWidth', 1.2); hold on;
plot(continuousTable.K, continuousTable.P2_imag, 'LineWidth', 1.2);
plot(continuousTable.K, continuousTable.P3_imag, 'LineWidth', 1.2);
plot(continuousTable.K, continuousTable.P4_imag, 'LineWidth', 1.2);
grid on; xlabel('K'); ylabel('Imag Part'); title('Pole Imaginary Parts'); xline(2.36, 'k:');

subplot(4, 3, 9);
plot(continuousTable.K, continuousTable.num_complex_pairs, 'm-', 'LineWidth', 1.2); grid on; xlabel('K'); ylabel('Pairs'); title('Number of Complex Pairs'); xline(2.36, 'k:');

subplot(4, 3, 10);
plot(continuousTable.K, continuousTable.PO_error_percent, 'm-', 'LineWidth', 1.2); grid on; xlabel('K'); ylabel('Error (%)'); title('PO Error'); xline(2.36, 'k:');

subplot(4, 3, 11);
plot(continuousTable.K, continuousTable.tp_error_percent, 'm-', 'LineWidth', 1.2); grid on; xlabel('K'); ylabel('Error (%)'); title('tp Error'); xline(2.36, 'k:');

subplot(4, 3, 12);
plot(continuousTable.K, continuousTable.ts_error_percent, 'm-', 'LineWidth', 1.2); grid on; xlabel('K'); ylabel('Error (%)'); title('ts Error'); xline(2.36, 'k:');

sgtitle('Continuous Gain Analysis Summary');
saveas(fig, fullfile(outputDir, 'fig_gain_analysis_comprehensive.png'));
close(fig);

% -------------------------------------------------------------------------
% DISCUSSION NOTES AS COMMENTS
% -------------------------------------------------------------------------
% WHY THE 2ND-ORDER APPROXIMATION FAILS FOR LOW K (K = 0.1 TO 2.36):
% 1. The closed-loop system still has two complex conjugate pairs, so there
%    is no clean single dominant mode.
% 2. The non-dominant oscillatory pair is not much faster, so it affects the
%    transient shape instead of dying out quickly.
% 3. The open-loop zeros at s = -1 and s = -0.2 act like derivative shaping
%    terms and strongly change overshoot and rise behavior.
% 4. The double integrator at the origin makes the dynamics fundamentally
%    different from the textbook second-order form.
% 5. The interaction of the two oscillatory modes suppresses the overshoot
%    that a single-pair model would usually predict.
%
% THE CRITICAL TRANSITION AT K APPROX 2.36:
% 1. Below this value there are two complex pairs.
% 2. Above this value there is one complex pair and two real poles.
% 3. This is where the root-locus branches from the complex open-loop poles
%    meet the real axis and split into real-axis branches.
% 4. The approximation error spikes because the dominant-pair definition is
%    changing at the same time as the pole geometry changes.
%
% WHY APPROXIMATION IMPROVES FOR HIGH K (K = 10 TO 50):
% 1. The system settles into one dominant complex pair plus two real poles.
% 2. The real poles move farther left, so they decay faster and become less
%    important in the transient response.
% 3. The complex pair becomes a meaningful reduced-order model.
% 4. The zeros still matter, so the approximation improves but never becomes
%    perfect.
%
% EFFECT OF POLE VARIATION:
% 1. A +/-20 percent change in the open-loop pole pair causes large changes
%    in overshoot and timing.
% 2. This demonstrates strong sensitivity to parametric uncertainty.
% 3. Robustness analysis is therefore essential.
%
% STEADY-STATE ERROR:
% 1. The system is Type 2 because of the double integrator at the origin.
% 2. For a unit step input, the theoretical steady-state error is zero.
% 3. The simulated value should therefore be very small.

disp('Experiment 1 completed. All CSV files and PNG figures have been saved.');

% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function rows = analyzeGainSweep(Klist, variantLabel)
    rows = repmat(struct(), numel(Klist), 1);
    for idx = 1:numel(Klist)
        K = Klist(idx);
        [polesSorted, dominantPole, zeta, omegaN, omegaD, numComplexPairs, numRealPoles] = analyzeClosedLoop(K, variantLabel);
        [~, ~, metrics] = simulateAndMeasureStep(K, variantLabel, polesSorted, dominantPole);

        rows(idx).K = K;
        rows(idx).zeta = zeta;
        rows(idx).omega_n = omegaN;
        rows(idx).omega_d = omegaD;
        rows(idx).num_complex_pairs = numComplexPairs;
        rows(idx).num_real_poles = numRealPoles;
        rows(idx).PO_measured = metrics.PO;
        rows(idx).PO_approx = metrics.PO_approx;
        rows(idx).PO_error_percent = safePercentError(metrics.PO, metrics.PO_approx);
        rows(idx).tp_measured = metrics.tp;
        rows(idx).tp_approx = metrics.tp_approx;
        rows(idx).tp_error_percent = safePercentError(metrics.tp, metrics.tp_approx);
        rows(idx).ts_measured = metrics.ts;
        rows(idx).ts_approx = metrics.ts_approx;
        rows(idx).ts_error_percent = safePercentError(metrics.ts, metrics.ts_approx);
        rows(idx).ess_measured = metrics.ess;
        rows(idx).y_ss = metrics.yss;
        rows(idx).P1_real = real(polesSorted(1));
        rows(idx).P1_imag = imag(polesSorted(1));
        rows(idx).P2_real = real(polesSorted(2));
        rows(idx).P2_imag = imag(polesSorted(2));
        rows(idx).P3_real = real(polesSorted(3));
        rows(idx).P3_imag = imag(polesSorted(3));
        rows(idx).P4_real = real(polesSorted(4));
        rows(idx).P4_imag = imag(polesSorted(4));
    end
end

function [polesSorted, dominantPole, zeta, omegaN, omegaD, numComplexPairs, numRealPoles] = analyzeClosedLoop(K, variantLabel)
    tol = 1e-7;
    clDen = closedLoopDenominator(K, variantLabel);
    poles = roots(clDen);
    polesSorted = sortPolesForReporting(poles);

    numRealPoles = sum(abs(imag(poles)) <= tol);
    numComplexPairs = sum(imag(poles) > tol);

    complexRoots = poles(imag(poles) > tol);
    if isempty(complexRoots)
        dominantPole = polesSorted(1);
        zeta = NaN;
        omegaN = NaN;
        omegaD = NaN;
        return;
    end

    [~, bestIdx] = max(real(complexRoots));
    dominantPole = complexRoots(bestIdx);
    sigma = -real(dominantPole);
    omegaD = abs(imag(dominantPole));
    omegaN = sqrt(sigma^2 + omegaD^2);
    if omegaN > 0
        zeta = sigma / omegaN;
    else
        zeta = NaN;
    end
end

function [t, y, metrics] = simulateAndMeasureStep(K, variantLabel, polesSorted, dominantPole)
    sys = closedLoopTF(K, variantLabel);
    tfinal = adaptiveFinalTime(polesSorted, dominantPole);
    nPts = max(800, min(2500, ceil(50 * tfinal)));
    t = linspace(0, tfinal, nPts);
    y = step(sys, t);
    y = squeeze(y);
    metrics = measureResponseMetrics(t, y, dominantPole);
end

function metrics = measureResponseMetrics(t, y, dominantPole)
    n = numel(y);
    tailStart = max(1, floor(0.9 * n));
    yss = mean(y(tailStart:end));
    if ~isfinite(yss) || abs(yss) < 1e-10
        yss = y(end);
    end
    if ~isfinite(yss) || abs(yss) < 1e-10
        yss = 1;
    end

    [peakVal, peakIdx] = max(y);
    tp = t(peakIdx);
    PO = max(0, (peakVal - yss) / abs(yss) * 100);

    band = 0.02 * abs(yss);
    outOfBand = find(abs(y - yss) > band);
    if isempty(outOfBand)
        ts = 0;
    else
        lastViolation = outOfBand(end);
        if lastViolation < n
            ts = t(lastViolation + 1);
        else
            ts = t(end);
        end
    end

    ess = abs(1 - yss);

    sigma = -real(dominantPole);
    omegaD = abs(imag(dominantPole));
    omegaN = sqrt(sigma^2 + omegaD^2);
    if isfinite(sigma) && isfinite(omegaN) && sigma > 0 && omegaN > 0 && sigma < omegaN
        zeta = sigma / omegaN;
        PO_approx = 100 * exp((-pi * zeta) / sqrt(max(1e-12, 1 - zeta^2)));
        tp_approx = pi / max(omegaD, eps);
        ts_approx = 4 / (zeta * omegaN);
    else
        PO_approx = NaN;
        tp_approx = NaN;
        ts_approx = NaN;
    end

    metrics = struct();
    metrics.PO = PO;
    metrics.tp = tp;
    metrics.ts = ts;
    metrics.ess = ess;
    metrics.yss = yss;
    metrics.PO_approx = PO_approx;
    metrics.tp_approx = tp_approx;
    metrics.ts_approx = ts_approx;
    metrics.ess_approx = 0;
end

function sys = closedLoopTF(K, variantLabel)
    den = closedLoopDenominator(K, variantLabel);
    num = K * [5 6 1];
    sys = tf(num, den);
end

function den = closedLoopDenominator(K, variantLabel)
    switch lower(string(variantLabel))
        case "nominal"
            a = 2;
            b = 10;
        case "negative20"
            a = 1.6;
            b = 6.4;
        case "positive20"
            a = 2.4;
            b = 12.96;
        otherwise
            error('Unknown variant label: %s', variantLabel);
    end
    den = [1, a, b + 5*K, 6*K, K];
end

function tfinal = adaptiveFinalTime(polesSorted, dominantPole)
    stablePoles = polesSorted(real(polesSorted) < 0);
    if isempty(stablePoles)
        stablePoles = polesSorted;
    end
    slowRates = abs(real(stablePoles(real(stablePoles) < 0)));
    if isempty(slowRates)
        slowRates = 0.1;
    end
    slowReal = min(slowRates);

    sigma = abs(real(dominantPole));
    if ~isfinite(sigma) || sigma <= 0
        sigma = slowReal;
    end

    omegaD = abs(imag(dominantPole));
    settleTime = 10 / sigma;
    if omegaD > 0
        periodTime = 6 * 2 * pi / omegaD;
    else
        periodTime = 0;
    end
    decayTime = 8 / slowReal;

    tfinal = max([15, settleTime, periodTime, decayTime]);
    tfinal = min(tfinal, 150);
end

function polesSorted = sortPolesForReporting(poles)
    tol = 1e-7;
    positiveImag = poles(imag(poles) > tol);
    realPoles = poles(abs(imag(poles)) <= tol);

    [~, idxPos] = sort(real(positiveImag), 'descend');
    positiveImag = positiveImag(idxPos);
    negativeImag = conj(positiveImag);

    [~, idxReal] = sort(real(realPoles), 'descend');
    realPoles = realPoles(idxReal);

    polesSorted = [positiveImag(:); negativeImag(:); realPoles(:)];
    if numel(polesSorted) < 4
        polesSorted = [polesSorted; complex(NaN, NaN) * ones(4 - numel(polesSorted), 1)];
    end
    polesSorted = polesSorted(1:4);
end

function err = safePercentError(measured, approx)
    if ~isfinite(measured) || ~isfinite(approx)
        err = NaN;
    elseif abs(measured) < 1e-12
        if abs(approx) < 1e-12
            err = 0;
        else
            err = NaN;
        end
    else
        err = abs(measured - approx) / abs(measured) * 100;
    end
end

function out = sanitizeValue(val)
    out = strrep(sprintf('%.6f', val), '.', '_');
end

function grid = buildCommonTimeGrid(timeVectors)
    maxT = 0;
    for i = 1:numel(timeVectors)
        maxT = max(maxT, timeVectors{i}(end));
    end
    nPts = max(1200, min(4000, ceil(80 * maxT)));
    grid = linspace(0, maxT, nPts);
end