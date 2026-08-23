%% ========================================================================
%  PART A -- ROOT-LOCUS CASE STUDIES & BREAKAWAY/BREAK-IN ANALYSIS
%  Reproduces report Figures 2-7 and Tables 1-2.
%
%  System:   G0(s) = (s+1) / [s^2(s+a)],   G(s) = K*G0(s),   H(s) = 1
%  Char. eq: s^3 + a*s^2 + K*s + K = 0
%
%  Requires: Control System Toolbox (tf, rlocus)
%  Needs on path: charEqCoeffs.m, viridisLike.m, breakawayPoints.m
% ========================================================================
clc; clear; close all;

outdir = 'figures_partA';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ---- Table 1: open-loop pole / asymptote centroid vs a --------------
a_cases = [1000 12 10 9 5 3 2 1];
fprintf('\n=== Table 1: Open-loop pole / asymptote centroid ===\n');
fprintf('%8s %14s %16s\n', 'a', 'third pole', 'centroid sigma_A');
for a = a_cases
    fprintf('%8g %14g %16.2f\n', a, -a, (1-a)/2);
end

%% ---- Figure 2: Root-locus overlay for several values of a -----------
a_moderate = [1 2 3 5 7 9 10 12 20];
cmap = viridisLike(numel(a_moderate));

figure('visible', 'off', 'Position', [50 50 850 650]);
hold on; box on; grid on;
legend_handles = [];
legend_labels  = {};
for i = 1:numel(a_moderate)
    a = a_moderate(i);
    G0 = tf([1 1], [1 a 0 0]);          % G0(s) = (s+1)/(s^2(s+a))
    [r, K] = rlocus(G0);                 % r: 3 x N pole trajectories
    h = plot(real(r).', imag(r).', 'Color', cmap(i,:), 'LineWidth', 1.6);
    legend_handles(end+1) = h(1); %#ok<SAGROW>
    legend_labels{end+1}  = sprintf('a = %g', a); %#ok<SAGROW>
    plot(-a, 0, 'x', 'Color', cmap(i,:), 'MarkerSize', 9, 'LineWidth', 2);
end
plot(-1, 0, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
xline(0, 'k-'); yline(0, 'k-');
xlim([-22 3]); ylim([-16 16]);
xlabel('Real Axis (s^{-1})'); ylabel('Imaginary Axis (s^{-1})');
title('Root Locus Overlay for Several Values of a  (K: 0 \rightarrow \infty)');
legend(legend_handles, legend_labels, 'Location', 'northwest', 'NumColumns', 2, 'FontSize', 8);
print(fullfile(outdir, 'fig02_rootlocus_overlay.png'), '-dpng', '-r150');

%% ---- Figure 3: Zoom near the critical transition at a = 9 -----------
a_zoom = [7 8 9 10 12];
figure('visible', 'off', 'Position', [50 50 1500 350]);
for i = 1:numel(a_zoom)
    a = a_zoom(i);
    subplot(1, numel(a_zoom), i);
    G0 = tf([1 1], [1 a 0 0]);
    [r, K] = rlocus(G0);
    plot(real(r).', imag(r).', 'Color', [0.17 0.36 0.54], 'LineWidth', 1.4); hold on;
    plot(-1, 0, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
    plot(-a, 0, 'x', 'Color', [0.76 0.27 0.05], 'MarkerSize', 7, 'LineWidth', 2);
    xline(0, 'k-'); yline(0, 'k-');
    xlim([-a-2 2]); ylim([-6 6]);
    if a < 9, tag = 'no loop'; elseif a == 9, tag = 'borderline'; else, tag = 'loop present'; end
    title(sprintf('a = %g\n(%s)', a, tag), 'FontSize', 10);
    xlabel('Re(s)'); grid on; box on;
    if i == 1, ylabel('Im(s)'); end
end
sgtitle('Locus Detail Near the Critical Transition at a = 9');
print(fullfile(outdir, 'fig03_rootlocus_zoom_a9.png'), '-dpng', '-r150');

%% ---- Figure 4: Root locus for a = 1000 (large-a extreme) ------------
a = 1000;
G0 = tf([1 1], [1 a 0 0]);
figure('visible', 'off', 'Position', [50 50 700 560]);
[r, K] = rlocus(G0);
h1 = plot(real(r).', imag(r).', 'Color', [0.17 0.36 0.54], 'LineWidth', 1.6); hold on;
h2 = plot(-1, 0, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
h3 = plot(-a, 0, 'x', 'Color', [0.76 0.27 0.05], 'MarkerSize', 9, 'LineWidth', 2);
h4 = plot(0, 0, 'kx', 'MarkerSize', 9, 'LineWidth', 2);
h5 = xline((1-a)/2, '--', 'Color', [0.2 0.6 0.2], 'LineWidth', 1);
xline(0, 'k-'); yline(0, 'k-');
xlim([-1150 150]); ylim([-320 320]);
xlabel('Real Axis (s^{-1})'); ylabel('Imaginary Axis (s^{-1})');
title(sprintf('Root Locus for a = %d (large-a extreme)', a));
legend([h1(1) h2 h3 h4 h5], {'locus', 'zero s=-1', sprintf('pole s=-%d', a), ...
    'double pole s=0', sprintf('asymptote Re(s)=%.1f', (1-a)/2)}, 'Location', 'east', 'FontSize', 9);
grid on; box on;
print(fullfile(outdir, 'fig04_rootlocus_a1000.png'), '-dpng', '-r150');

%% ---- Figure 5: Boundary case a = 1 (pole-zero cancellation) ---------
a = 1;
G0 = tf([1 1], [1 a 0 0]);   % pole/zero at s=-1 cancel algebraically in G0
figure('visible', 'off', 'Position', [50 50 620 540]);
[r, K] = rlocus(G0);
h1 = plot(real(r).', imag(r).', 'Color', [0.17 0.36 0.54], 'LineWidth', 1.6); hold on;
h2 = plot(-1, 0, 'o', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'MarkerSize', 9);
h3 = plot(0, 0, 'kx', 'MarkerSize', 9, 'LineWidth', 2);
xline(0, 'k-'); yline(0, 'k-');
xlim([-1.3 0.3]); ylim([-6.5 6.5]);
xlabel('Real Axis (s^{-1})'); ylabel('Imaginary Axis (s^{-1})');
title({'Root Locus for the Boundary Case a = 1', 'reduces to G(s)=K/s^2, locus on j\omega-axis'});
legend([h2 h3], {'coincident pole/zero (s=-1)', 'double pole s=0'}, 'Location', 'southeast', 'FontSize', 9);
grid on; box on;
print(fullfile(outdir, 'fig05_rootlocus_a1_boundary.png'), '-dpng', '-r150');

%% ---- Figure 6: Instability demonstration, 0 < a < 1 ------------------
a = 0.5;
G0 = tf([1 1], [1 a 0 0]);
figure('visible', 'off', 'Position', [50 50 700 560]);
hp = patch([0 6 6 0], [-4 -4 4 4], 'r', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
hold on;
[r, K] = rlocus(G0);
h1 = plot(real(r).', imag(r).', 'Color', [0.17 0.36 0.54], 'LineWidth', 1.6);
h2 = plot(-1, 0, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
h3 = plot(-a, 0, 'x', 'Color', [0.76 0.27 0.05], 'MarkerSize', 9, 'LineWidth', 2);
h4 = plot(0, 0, 'kx', 'MarkerSize', 9, 'LineWidth', 2);
xline(0, 'k-', 'LineWidth', 1.1); yline(0, 'k-');
xlim([-3.5 3.5]); ylim([-4 4]);
xlabel('Real Axis (s^{-1})'); ylabel('Imaginary Axis (s^{-1})');
title(sprintf('Root Locus for a = %.1f (0<a<1) -- branches enter RHP for all K>0', a));
legend([hp h2 h3 h4], {'Right-half plane (unstable)', 'zero s=-1', sprintf('pole s=-%.1f', a), 'double pole s=0'}, ...
    'Location', 'northwest', 'FontSize', 8.5);
grid on; box on;
print(fullfile(outdir, 'fig06_rootlocus_instability.png'), '-dpng', '-r150');

%% ---- Figure 7: Breakaway/break-in discriminant D(a)=(a-1)(a-9) ------
a_line = linspace(0.2, 15, 600);
D_line = (a_line + 3).^2 - 16*a_line;

figure('visible', 'off', 'Position', [50 50 750 480]);
plot(a_line, D_line, 'LineWidth', 2.2, 'Color', [0.17 0.36 0.54]); hold on;
yline(0, 'k-');
xline(1, '--', 'Color', [0.76 0.27 0.05]);
xline(9, '--', 'Color', [0.76 0.27 0.05]);
fill_pos = D_line; fill_pos(fill_pos < 0) = NaN;
fill_neg = D_line; fill_neg(fill_neg > 0) = NaN;
area(a_line, fill_pos, 'FaceColor', [0.56 0.75 0.42], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
area(a_line, fill_neg, 'FaceColor', [0.88 0.56 0.43], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
xlabel('Open-loop pole parameter, a');
ylabel('Discriminant  D(a) = (a+3)^2 - 16a = (a-1)(a-9)');
title('Breakaway/Break-in Discriminant vs. a -- origin of the critical value a=9');
text(1.6, 45, 'a=1 (pole-zero cancellation)', 'FontSize', 8.5);
text(9.6, 45, 'a=9 (loop shrinks to a point)', 'FontSize', 8.5);
grid on; box on;
print(fullfile(outdir, 'fig07_discriminant.png'), '-dpng', '-r150');

%% ---- Table 2: breakaway/break-in points for a > 9 --------------------
fprintf('\n=== Table 2: Breakaway/break-in points (a > 9) ===\n');
fprintf('%6s %10s %14s %14s %14s %14s\n', 'a', 'D(a)', 'break-in s', 'K(break-in)', 'breakaway s', 'K(breakaway)');
for a = [9 10 12 20 50]
    [D, s_bi, K_bi, s_ba, K_ba] = breakawayPoints(a);
    if a == 9
        fprintf('%6g %10.1f %14s %14.2f %14s %14.2f\n', a, D, '(coincident,s=-3)', K_bi, '(coincident,s=-3)', K_ba);
    else
        fprintf('%6g %10.1f %14.3f %14.2f %14.3f %14.2f\n', a, D, s_bi, K_bi, s_ba, K_ba);
    end
end

fprintf('\nPart A complete. Figures saved to %s/\n', outdir);
