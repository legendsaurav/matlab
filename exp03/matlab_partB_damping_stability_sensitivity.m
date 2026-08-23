%% ========================================================================
%  PART B -- DAMPING-RATIO APPROXIMATION, ROUTH-HURWITZ, SENSITIVITY
%  Reproduces report Figures 8, 9, 21 and Tables 3, 4, 5, 10.
%
%  System:   s^3 + a*s^2 + K*s + K = 0
%  Approx:   zeta ~= (a-1)/(2*sqrt(K))   for K >> (a-1)^2/4
%
%  Requires: base MATLAB only (roots, polyfit) -- no toolboxes needed.
%  Needs on path: charEqCoeffs.m, dominantComplexPole.m
% ========================================================================
clc; clear; close all;

outdir = 'figures_partB';
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ---- Exact damping ratio table: a = 1..50, five gains ----------------
a_values = (1:50)';
K_values = [10000 20000 30000 50000 70000];

zeta_exact  = zeros(numel(a_values), numel(K_values));
zeta_approx = zeros(size(zeta_exact));
for j = 1:numel(K_values)
    K = K_values(j);
    for i = 1:numel(a_values)
        a = a_values(i);
        s = dominantComplexPole(a, K);
        zeta_exact(i,j)  = -real(s) / abs(s);
        zeta_approx(i,j) = (a - 1) / (2*sqrt(K));
    end
end

fprintf('\n=== Damping ratio (a=1..9 excerpt, matches report Table 2) ===\n');
fprintf('%4s', 'a'); fprintf('%14d', K_values); fprintf('\n');
for i = 1:9
    fprintf('%4d', a_values(i));
    fprintf('%14.6f', zeta_exact(i,:));
    fprintf('\n');
end

%% ---- Figure 8: exact vs approximate damping ratio --------------------
cmap = plasmaLike(numel(K_values));
figure('visible', 'off', 'Position', [50 50 760 540]);
hold on; grid on; box on;
hE = gobjects(numel(K_values),1);
for j = 1:numel(K_values)
    hE(j) = plot(a_values, zeta_exact(:,j), 'Color', cmap(j,:), 'LineWidth', 2.1);
    plot(a_values, zeta_approx(:,j), '--', 'Color', cmap(j,:), 'LineWidth', 1.3);
end
xlabel('Open-Loop Pole Position, a'); ylabel('Damping Ratio, \zeta');
title('Damping Ratio vs. a -- Exact (roots) vs. Approximate Formula');
legend(hE, arrayfun(@(K) sprintf('K = %d (exact)', K), K_values, 'UniformOutput', false), ...
    'Location', 'northwest', 'FontSize', 8.5);
print(fullfile(outdir, 'fig08_damping_exact_vs_approx.png'), '-dpng', '-r150');

%% ---- Figure 9: approximation residual --------------------------------
figure('visible', 'off', 'Position', [50 50 740 500]);
hold on; grid on; box on;
for j = 1:numel(K_values)
    plot(a_values, zeta_exact(:,j) - zeta_approx(:,j), 'Color', cmap(j,:), 'LineWidth', 1.9);
end
yline(0, 'k-');
xlabel('Open-Loop Pole Position, a'); ylabel('Residual  \zeta_{exact} - \zeta_{approx}');
title('Approximation Error vs. a (shrinks as K grows, as predicted)');
legend(arrayfun(@(K) sprintf('K = %d', K), K_values, 'UniformOutput', false), ...
    'Location', 'northwest', 'FontSize', 9);
print(fullfile(outdir, 'fig09_damping_approx_error.png'), '-dpng', '-r150');

%% ---- Table 3: approximation error summary ----------------------------
% NOTE: at a=1, zeta_exact is (numerically) ~0, so the raw percentage
% error 100*|err|/zeta_exact blows up / divides by ~0. We exclude that
% single degenerate point from the percentage-error mean, exactly as the
% underlying analysis does -- this matters: computing it unguarded
% contaminates the whole mean with NaN/Inf.
fprintf('\n=== Table 3: Approximation error summary ===\n');
fprintf('%8s %12s %12s %12s\n', 'K', 'RMSE', 'MaxAbsErr', 'Mean pct Err');
for j = 1:numel(K_values)
    err = zeta_exact(:,j) - zeta_approx(:,j);
    rmse = sqrt(mean(err.^2));
    maxerr = max(abs(err));
    valid = zeta_exact(:,j) > 1e-9;         % excludes the a=1 degenerate point
    pcterr = mean(100*abs(err(valid)) ./ zeta_exact(valid,j));
    fprintf('%8d %12.6f %12.6f %11.3f%%\n', K_values(j), rmse, maxerr, pcterr);
end

%% ---- Table 4: linear fit vs theoretical slope/intercept --------------
fprintf('\n=== Table 4: Linear fit zeta = m*a + c  vs. theory ===\n');
fprintf('%8s %12s %14s %14s %16s %10s\n', 'K', 'fit slope', 'theory 1/2sqrtK', 'fit intercept', 'theory -1/2sqrtK', 'R^2');
for j = 1:numel(K_values)
    K = K_values(j);
    p = polyfit(a_values, zeta_exact(:,j), 1);
    m_theory = 1/(2*sqrt(K));
    yfit = polyval(p, a_values);
    ss_res = sum((zeta_exact(:,j) - yfit).^2);
    ss_tot = sum((zeta_exact(:,j) - mean(zeta_exact(:,j))).^2);
    R2 = 1 - ss_res/ss_tot;
    fprintf('%8d %12.6f %14.6f %14.6f %16.6f %10.6f\n', K, p(1), m_theory, p(2), -m_theory, R2);
end

%% ---- Routh-Hurwitz numerical verification (Table 5) ------------------
% Routh array for s^3+a*s^2+K*s+K=0:
%   s^3 |   1        K
%   s^2 |   a        K
%   s^1 | K(a-1)/a   0
%   s^0 |   K
fprintf('\n=== Table 5: Routh-Hurwitz verification (K = 10000) ===\n');
fprintf('%8s %16s %14s\n', 'a', 's^1 row K(a-1)/a', 'numeric check');
test_a = [0.2 0.5 0.8 0.95 1 1.05 1.5 2 5 10 50];
for a = test_a
    K = 10000;
    s1_entry = K*(a-1)/a;
    r = roots(charEqCoeffs(a, K));
    if all(real(r) < -1e-9)
        verdict = 'stable';
    elseif any(real(r) > 1e-9)
        verdict = 'unstable';
    else
        verdict = 'marginal';
    end
    fprintf('%8.2f %16.3f %14s\n', a, s1_entry, verdict);
end
% Repeat across K to confirm the verdict is K-independent for each a
fprintf('\nCross-check across K = [1 10 100 1000 10000] (same verdict every time):\n');
all_match = true;
Kset = [1 10 100 1000 10000];
for a = test_a
    verdicts = cell(1, numel(Kset));    % cell array of char vectors: works
    for kk = 1:numel(Kset)              % identically on every MATLAB/Octave
        r = roots(charEqCoeffs(a, Kset(kk)));
        if all(real(r) < -1e-9), verdicts{kk} = 'stable';
        elseif any(real(r) > 1e-9), verdicts{kk} = 'unstable';
        else, verdicts{kk} = 'marginal'; end
    end
    if numel(unique(verdicts)) > 1
        all_match = false;
        fprintf('  a=%.2f: MISMATCH across K -> %s\n', a, strjoin(verdicts, ','));
    end
end
if all_match
    fprintf('  Confirmed: verdict depends only on a, never on K, for every a tested.\n');
end

%% ---- Pole sensitivity ds/da: analytical vs finite-difference ---------
K_fixed_list = [10000 30000 70000];
da = 1e-4;
a_sens = 2:2:50;

fprintf('\n=== Table 10: Sensitivity ds/da, analytical vs finite-difference (K=10000) ===\n');
fprintf('%6s %24s %24s %12s\n', 'a', 'analytical', 'finite-diff', '|diff|');

sens_analytic = zeros(numel(a_sens), numel(K_fixed_list));
sens_numeric  = zeros(size(sens_analytic));
for j = 1:numel(K_fixed_list)
    K = K_fixed_list(j);
    for i = 1:numel(a_sens)
        a = a_sens(i);
        s0 = dominantComplexPole(a, K);
        s_analytic = -(s0^2) / (3*s0^2 + 2*a*s0 + K);

        rp = roots(charEqCoeffs(a+da, K));
        [~, ip] = min(abs(rp - s0));
        rm = roots(charEqCoeffs(a-da, K));
        [~, im] = min(abs(rm - s0));
        s_numeric = (rp(ip) - rm(im)) / (2*da);

        sens_analytic(i,j) = abs(s_analytic);
        sens_numeric(i,j)  = abs(s_numeric);

        if j == 1 && ismember(a, [2 10 20 30 50])
            fprintf('%6d %11.5f%+.5fi %11.5f%+.5fi %12.2e\n', a, ...
                real(s_analytic), imag(s_analytic), real(s_numeric), imag(s_numeric), ...
                abs(s_analytic - s_numeric));
        end
    end
end

%% ---- Figure 18: pole sensitivity vs a ---------------------------------
figure('visible', 'off', 'Position', [50 50 720 520]);
hold on; grid on; box on;
cmap2 = viridisLike(numel(K_fixed_list));
for j = 1:numel(K_fixed_list)
    plot(a_sens, sens_analytic(:,j), 'o-', 'Color', cmap2(j,:), 'MarkerSize', 4, 'LineWidth', 1.6);
end
xlabel('Open-Loop Pole Position, a'); ylabel('Pole Sensitivity  |\partial s/\partial a|');
title('Closed-Loop Pole Sensitivity to the Parameter a');
legend(arrayfun(@(K) sprintf('K = %d', K), K_fixed_list, 'UniformOutput', false), 'Location', 'northwest', 'FontSize', 9.5);
print(fullfile(outdir, 'fig21_pole_sensitivity.png'), '-dpng', '-r150');

fprintf('\nPart B complete. Figures saved to %s/\n', outdir);
