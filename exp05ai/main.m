function results = main()
%% MAIN  Canonical analysis entry point for the continuous-time model.
% Model from the handwritten specification:
%   G(s) = 10 (s + 2.5)^2 / ((s + 10) (s^2 + 0.12))
%
% All generated figures are written to assets/ and all tabular data to
% artifacts/. Run this file from any current folder with: main

clc;
close all;

projectRoot = fileparts(mfilename('fullpath'));
assetDir = fullfile(projectRoot, 'assets');
artifactDir = fullfile(projectRoot, 'artifacts');
if ~exist(assetDir, 'dir'), mkdir(assetDir); end
if ~exist(artifactDir, 'dir'), mkdir(artifactDir); end

s = tf('s');
G = 10 * (s + 2.5)^2 / ((s + 12) * (s^2 + 0.12));

fprintf('Continuous-time model\n');
disp(G);

% ---- Pole-zero data ---------------------------------------------------
plantPoles = pole(G);
plantZeros = zero(G);
T = table([plantZeros; plantPoles], ...
    [repmat("zero", numel(plantZeros), 1); repmat("pole", numel(plantPoles), 1)], ...
    'VariableNames', {'value', 'type'});
writetable(T, fullfile(artifactDir, 'table1_poles_zeros.csv'));

fig = figure('Visible', 'off', 'Position', [100 100 720 560]);
pzmap(G);
grid on;
title('Pole-Zero Map of G(s)');
exportgraphics(fig, fullfile(assetDir, 'fig1_polezero.png'), 'Resolution', 220);
close(fig);

% ---- Open-loop frequency response ------------------------------------
frequency = logspace(-3, 3, 3000);
[mag, phase] = bode(G, frequency);
mag = squeeze(mag);
phase = squeeze(phase);
T = table(frequency(:), 20 * log10(mag(:)), phase(:), ...
    'VariableNames', {'frequency_rad_s', 'magnitude_dB', 'phase_deg'});
writetable(T, fullfile(artifactDir, 'table15_bode_plant.csv'));

fig = figure('Visible', 'off', 'Position', [100 100 820 650]);
bodeplot(G, frequency);
grid on;
title('Bode Plot of G(s)');
exportgraphics(fig, fullfile(assetDir, 'fig20_bode_plant.png'), 'Resolution', 220);
close(fig);

% ---- Closed-loop gain sweep -----------------------------------------
gain = logspace(-3, 3, 500);
closedLoopPoles = zeros(numel(gain), numel(plantPoles));
stable = false(size(gain));
for index = 1:numel(gain)
    polesAtGain = pole(feedback(gain(index) * G, 1));
    closedLoopPoles(index, :) = polesAtGain(:).';
    stable(index) = all(real(polesAtGain) < 0);
end
T = table(gain(:), stable(:), ...
    'VariableNames', {'K', 'stable'});
for index = 1:size(closedLoopPoles, 2)
    T.(sprintf('pole%d_real', index)) = real(closedLoopPoles(:, index));
    T.(sprintf('pole%d_imag', index)) = imag(closedLoopPoles(:, index));
end
writetable(T, fullfile(artifactDir, 'table2_stability_sweep.csv'));

fig = figure('Visible', 'off', 'Position', [100 100 850 600]);
rlocus(G);
grid on;
title('Root Locus of G(s)');
exportgraphics(fig, fullfile(assetDir, 'fig2_rootlocus_stability.png'), 'Resolution', 220);
close(fig);

% ---- Representative closed-loop response and margins ----------------
Kdesign = 1;
Tclosed = feedback(Kdesign * G, 1);
[stepResponse, stepTime] = step(Tclosed);
writetable(table(stepTime(:), stepResponse(:), ...
    'VariableNames', {'time_s', 'response'}), ...
    fullfile(artifactDir, 'table_waterfall_stepresponses.csv'));

fig = figure('Visible', 'off', 'Position', [100 100 820 520]);
step(Tclosed);
grid on;
title(sprintf('Unity-Feedback Step Response, K = %.3g', Kdesign));
exportgraphics(fig, fullfile(assetDir, 'fig10_waterfall.png'), 'Resolution', 220);
close(fig);

[gainMargin, phaseMargin, ~, ~] = margin(G);
writetable(table(Kdesign, gainMargin, phaseMargin, ...
    'VariableNames', {'K_design', 'gain_margin', 'phase_margin_deg'}), ...
    fullfile(artifactDir, 'table17_margins_summary.csv'));

results = struct('root', projectRoot, 'assets', assetDir, ...
    'artifacts', artifactDir, 'G', G, 'closedLoop', Tclosed, ...
    'KDesign', Kdesign, 'stableGainCount', nnz(stable));
fprintf('Analysis complete. Figures: %s\nData: %s\n', assetDir, artifactDir);
end
