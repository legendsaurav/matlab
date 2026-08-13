%% EE208 Control Systems Lab Experiment 1 - Tail Outputs
% This helper script generates the pole variation analysis and the
% comprehensive gain-analysis plot from the saved continuous CSV.

outputDir = fileparts(mfilename('fullpath'));
continuousTable = readtable(fullfile(outputDir, 'gain_analysis_continuous.csv'));

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

for vIdx = 1:numel(poleVariations)
    variation = poleVariations{vIdx};
    [polesSorted, dominantPole] = analyzeClosedLoopTail(5, variation.label);
    [t, y, metrics] = simulateAndMeasureStepTail(5, variation.label, polesSorted, dominantPole);

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

    fieldPrefix = sanitizeTextTail(variation.label);
    variationSeries.([fieldPrefix '_t']) = t;
    variationSeries.([fieldPrefix '_y']) = y;
    variationPoles.(variation.label) = polesSorted;
end

variationTable = cell2table(variationRows, 'VariableNames', {'Variation_Label', 'K', 'PO', 'tp', 'ts', 'ess', 'Pole1_Real', 'Pole1_Imag', 'Pole2_Real', 'Pole2_Imag', 'Pole3_Real', 'Pole3_Imag', 'Pole4_Real', 'Pole4_Imag'});
writetable(variationTable, fullfile(outputDir, 'pole_variation.csv'));

timeGrid = buildCommonTimeGridTail({variationSeries.Negative20_t, variationSeries.Nominal_t, variationSeries.Positive20_t});
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

disp('Tail analysis complete.');

% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function [polesSorted, dominantPole] = analyzeClosedLoopTail(K, variantLabel)
    poles = roots(closedLoopDenominatorTail(K, variantLabel));
    polesSorted = sortPolesForReportingTail(poles);
    complexRoots = poles(imag(poles) > 1e-7);
    if isempty(complexRoots)
        dominantPole = polesSorted(1);
    else
        [~, bestIdx] = max(real(complexRoots));
        dominantPole = complexRoots(bestIdx);
    end
end

function [t, y, metrics] = simulateAndMeasureStepTail(K, variantLabel, polesSorted, dominantPole)
    sys = tf(K * [5 6 1], closedLoopDenominatorTail(K, variantLabel));
    tfinal = adaptiveFinalTimeTail(polesSorted, dominantPole);
    nPts = max(800, min(2500, ceil(50 * tfinal)));
    t = linspace(0, tfinal, nPts);
    y = squeeze(step(sys, t));
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
        ts = t(min(lastViolation + 1, n));
    end
    metrics = struct('PO', PO, 'tp', tp, 'ts', ts, 'ess', abs(1 - yss));
end

function den = closedLoopDenominatorTail(K, variantLabel)
    switch lower(string(variantLabel))
        case "nominal"
            a = 2; b = 10;
        case "negative20"
            a = 1.6; b = 6.4;
        case "positive20"
            a = 2.4; b = 12.96;
        otherwise
            error('Unknown variant label: %s', variantLabel);
    end
    den = [1, a, b + 5*K, 6*K, K];
end

function tfinal = adaptiveFinalTimeTail(polesSorted, dominantPole)
    stablePoles = polesSorted(real(polesSorted) < 0);
    if isempty(stablePoles)
        stablePoles = polesSorted;
    end
    slowRates = abs(real(stablePoles(real(stablePoles) < 0)));
    if isempty(slowRates)
        slowRates = 0.1;
    end
    sigma = max(1e-6, abs(real(dominantPole)));
    omegaD = abs(imag(dominantPole));
    settleTime = 10 / sigma;
    periodTime = 0;
    if omegaD > 0
        periodTime = 6 * 2 * pi / omegaD;
    end
    decayTime = 8 / min(slowRates);
    tfinal = min(150, max([15, settleTime, periodTime, decayTime]));
end

function polesSorted = sortPolesForReportingTail(poles)
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

function out = sanitizeTextTail(txt)
    out = regexprep(char(string(txt)), '[^A-Za-z0-9]', '_');
end

function grid = buildCommonTimeGridTail(timeVectors)
    maxT = 0;
    for i = 1:numel(timeVectors)
        maxT = max(maxT, timeVectors{i}(end));
    end
    grid = linspace(0, maxT, max(1200, min(4000, ceil(80 * maxT))));
end