%% ========================================================================
%  PART A -- ROOT-LOCUS CASE STUDIES & BREAKAWAY/BREAK-IN ANALYSIS
%  Reproduces report Figures 2-7 and Tables 1-2.
%
%  System:   G0(s) = (s+2.5) / [s^2(s+a)],   G(s) = K*G0(s),   H(s) = 1
%  Char. eq: s^3 + a*s^2 + K*s + 2.5*K = 0
%
%  Requires: Control System Toolbox (tf, rlocus)
%  Needs on path: charEqCoeffs.m, viridisLike.m, breakawayPoints.m
% ========================================================================
clc; clear; close all;

baseDir = fileparts(mfilename('fullpath'));
outdir = fullfile(baseDir, 'assets', 'figures_partA');
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ---- Table 1: open-loop pole / asymptote centroid vs a --------------
a_cases = [1000 12 10 9 5 3 2 1];
fprintf('\n=== Table 1: Open-loop pole / asymptote centroid ===\n');
fprintf('%8s %14s %16s\n', 'a', 'third pole', 'centroid sigma_A');
for a = a_cases
    fprintf('%8g %14g %16.2f\n', a, -a, (1-a)/2);
end

%% ---- Figure 2: Root-locus overlay for several values of a -----------
a_moderate = [2 3 5 7 9 10 12 20,22,25,30,34,38];
cmap = viridisLike(numel(a_moderate));

figure('visible', 'off', 'Position', [50 50 850 650]);
hold on; box on; grid on;
legend_handles = [];
legend_labels  = {};
for i = 1:numel(a_moderate)
    a = a_moderate(i);
    G0 = tf([1 2.5], [1 a 0 0]);          % G0(s) = (s+2.5)/(s^2(s+a))
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
% Search for loop/closes near a in [22,22.5] and produce zoomed plots verifying closure
a_search = 2.5




:0.005:100;   % fine search in suspected interval
crit_a = [];
closed_flags = false(size(a_search));
for idx = 1:numel(a_search)
    a_try = a_search(idx);
    G0_try = tf([1 2.5], [1 a_try 0 0]);
    % evaluate rlocus densely
    [r_try, K_try] = rlocus(G0_try, logspace(-8,8,5001));
    % For each K column, check if any conjugate pair forms a closed loop:
    % criterion: two poles form a complex-conjugate pair that for increasing K
    % move from LHP to RHP and then return (i.e., trajectory encircles a point).
    % Simpler practical test: check for presence of a closed arc: examine
    % whether imaginary part vs real part of a branch reverses direction (sign change in dRe/dK)
    closed_detected = false;
    for col = 1:size(r_try,2)
        p = r_try(:,col);
        % look for complex pair (nonzero imag) and non-monotonic real part along locus
        imcount = sum(abs(imag(p))>1e-4);
        if imcount >= 2
            % sort by real part magnitude to follow branch roughly
            re_vals = real(p(abs(imag(p))>1e-4));
            % check if real parts along K have non-monotonic behaviour
            if max(re_vals)-min(re_vals) > 1e-3
                % compute discrete derivative along K for that branch approximation
                % here approximate by checking successive K columns (across columns),
                % track one representative pole index by nearest matching across columns
                % Build trajectories by matching nearest poles between successive K steps
                trajRe = [];
                trajIm = [];
                % initialize with poles at first column
                prev = r_try(:,1);
                for c = 1:size(r_try,2)
                    curr = r_try(:,c);
                    % match each prev pole to closest curr pole
                    [~, loc] = min(abs(curr - prev(1)));
                    trajRe(end+1) = real(curr(loc)); %#ok<SAGROW>
                    trajIm(end+1) = imag(curr(loc)); %#ok<SAGROW>
                    prev = curr;
                end
                dre = diff(trajRe);
                % if real part changes sign in derivative (non-monotonic) we infer loop-like behavior
                if any(dre(1:end-1).*dre(2:end) < 0)
                    closed_detected = true;
                    break;
                end
            end
        end
    end
    closed_flags(idx) = closed_detected;
    if closed_detected
        crit_a = a_try;
        break;
    end
end

if isempty(crit_a)
    warning('No closed loop detected in [22,22.5]. Using coarse candidates for visualization.');
    % pick representative values for visualization around 22.497..22.5
    a_zoom = [22.497 22.499 22.5];
else
    % build small neighborhood around detected critical a for plotting
    a_zoom = unique(max(0, crit_a + (-0.003:0.001:0.003)));
end

% Create verification plots: show locus and annotate whether loop detected
figure('visible', 'off', 'Position', [50 50 1500 350]);
for i = 1:numel(a_zoom)
    a = a_zoom(i);
    G0 = tf([1 2.5], [1 a 0 0]);
    [r, K] = rlocus(G0, logspace(-8,8,2001));
    subplot(1, numel(a_zoom), i);
    plot(real(r).', imag(r).', 'Color', [0.17 0.36 0.54], 'LineWidth', 1.4); hold on;
    plot(-1, 0, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
    plot(-a, 0, 'x', 'Color', [0.76 0.27 0.05], 'MarkerSize', 7, 'LineWidth', 2);
    xline(0, 'k-'); yline(0, 'k-');
    xlim([-a-2 2]); ylim([-6 6]);
    % determine closed-loop flag for this specific a (reuse earlier check if available)
    idx = find(abs(a_search - a) < 1e-12, 1);
    if isempty(idx)
        % perform local check for this a
        flag = false;
        for col = 1:size(r,2)
            p = r(:,col);
            if sum(abs(imag(p))>1e-4) >= 2
                % simple non-monotonic real part test along K for one representative trajectory
                prev = r(:,1);
                trajRe = zeros(1,size(r,2));
                for c = 1:size(r,2)
                    curr = r(:,c);
                    [~, loc] = min(abs(curr - prev(1)));
                    trajRe(c) = real(curr(loc));
                    prev = curr;
                end
                dre = diff(trajRe);
                if any(dre(1:end-1).*dre(2:end) < 0)
                    flag = true; break;
                end
            end
        end
    else
        flag = closed_flags(idx);
    end
    tag = ternary(flag, 'loop present', 'no loop');
    title(sprintf('a = %.6g\n(%s)', a, tag), 'FontSize', 10);
    xlabel('Re(s)'); grid on; box on;
    if i == 1, ylabel('Im(s)'); end
end
sgtitle('Verification: Locus Detail Near Suspected Transition (22 \le a \le 22.5)');

% helper: simple inline ternary
function out = ternary(cond, a_true, a_false) %#ok<DEFNU>
    if cond, out = a_true; else out = a_false; end
end

print(fullfile(outdir, 'fig03_rootlocus_zoom_a22to225.png'), '-dpng', '-r150');

%% ---- Figure 4: Root locus for a = 1000 (large-a extreme) ------------
a = 1000;
G0 = tf([1 2.5], [1 a 0 0]);
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
D_line = (a_line - 2.5).*(a_line - 22.5);

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
ylabel('Discriminant  D(a) = (a+7.5)^2 - 40a = (a-2.5)(a-22.5)');
title('Breakaway/Break-in Discriminant vs. a -- critical values at a=2.5 and a=22.5');
text(3.0, 45, 'a=2.5 (critical breakaway)', 'FontSize', 8.5);
text(22.5, 45, 'a=22.5 (critical break-in)', 'FontSize', 8.5);
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
