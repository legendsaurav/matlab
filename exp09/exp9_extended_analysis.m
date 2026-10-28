%% EE208 Experiment 9 -- Extended Analysis (Section 7 of the report)
% (a) Extended root-locus verification of closed-loop stability
% (b) Continuous-domain gain (K) -> percent-overshoot (PO) study,
%     built on the 2nd-order anti-aliasing filter designed in
%     exp9_core_analysis.m (wc = 1.1547 rad/s, zeta0 = 1/sqrt(2))
% Requires base MATLAB + Control System Toolbox.

clear; clc; close all;

%% 7.1 Extended root-locus verification (k in [-300, 300])
k_vals = linspace(-300, 300, 12001);
max_abs_root = zeros(size(k_vals));

for i = 1:length(k_vals)
    k = k_vals(i);
    r = roots([1  k  1.7*k  -(k+1)]);   % closed-loop poles for this k
    max_abs_root(i) = max(abs(r));      % worst-case pole magnitude
end

[best_val, idx] = min(max_abs_root);    % best achievable case over the sweep
best_k = k_vals(idx);
fprintf('min_k max|root(k)| = %.6f at k = %.4f\n', best_val, best_k);
fprintf('Any strictly stable k found? %d\n', any(max_abs_root < 1));

figure('Name','Extended Root Locus');
subplot(1,2,1);
plot(real(roots([1 0 0 -1])), imag(roots([1 0 0 -1])), 'kx', 'MarkerSize', 10); hold on;
th = linspace(0,2*pi,400); plot(cos(th), sin(th), 'k--'); grid on; axis equal;
title('Root locus near |k| \leq 15'); xlabel('Real'); ylabel('Imag');
subplot(1,2,2);
plot(k_vals, max_abs_root); grid on; yline(1, 'r--');
xlabel('k'); ylabel('max_i |z_i(k)|'); title('Largest closed-loop pole magnitude vs. k');

%% 7.3 Closed-loop model and analytical overshoot formula
wc    = 1.1547;              % rad/s  (from the core-analysis filter design, Ts = 1 s)
zeta0 = 1/sqrt(2);           % nominal open-loop damping ratio (Butterworth)

zeta_cl = @(K,f) (zeta0*(1+f/100)) ./ sqrt(1+K);     % closed-loop damping ratio
PO      = @(K,f) 100*exp(-pi*zeta_cl(K,f)./sqrt(1-zeta_cl(K,f).^2));  % percent overshoot

%% 7.4 Step-response waterfall (K = 0.15 -> 300) + simulated overshoot
K_list = logspace(log10(0.15), log10(300), 16);
t = linspace(0, 20, 2500);

PO_sim = zeros(size(K_list));
Y = zeros(length(K_list), length(t));
for i = 1:length(K_list)
    K = K_list(i);
    num_T = K*wc^2;
    den_T = [1, sqrt(2)*wc, wc^2*(1+K)];
    y = step(tf(num_T, den_T), t);       % simulated closed-loop step response
    Y(i,:) = y;
    PO_sim(i) = (max(y)-y(end))/y(end)*100;   % overshoot measured from the simulation
end

figure('Name','Step-Response Waterfall');
waterfall(t, log10(K_list), Y); colormap('viridis');
xlabel('Time (s)'); ylabel('log_{10}(K)'); zlabel('y(t)');
title('Continuous Step-Response Waterfall, K = 0.15 \rightarrow 300');

%% 7.5 Approximate gain -> overshoot relationship
K_curve = logspace(log10(0.15), log10(300), 400);
PO_curve = PO(K_curve, 0);

figure('Name','Gain vs Overshoot');
semilogx(K_curve, PO_curve, 'LineWidth', 1.5); hold on;
semilogx(K_list, PO_sim, 'o', 'MarkerSize', 6);
grid on; xlabel('K (log scale)'); ylabel('Percent overshoot, PO (%)');
legend('Analytical PO(K)', 'Simulated step response', 'Location', 'northwest');
title('Approximate Gain \rightarrow Overshoot Relationship');

max_err = max(abs(PO(K_list,0) - PO_sim));
fprintf('Max |analytical - simulated| PO error: %.4f pp\n', max_err);

%% 7.6 Sensitivity to open-loop damping perturbation
K_grid = logspace(log10(0.15), log10(300), 140);
f_grid = linspace(-30, 30, 121);
[KK, FF] = meshgrid(K_grid, f_grid);

PO_grid = PO(KK, FF);     % fully vectorised -- no explicit loop required

figure('Name','Overshoot Sensitivity');
pcolor(KK, FF, PO_grid); shading interp; set(gca,'XScale','log'); colorbar;
xlabel('K (log scale)'); ylabel('Open-loop damping-ratio perturbation, f (%)');
title('Overshoot Sensitivity PO(K,f)');

% Table 4 anchor values (K = 1,5,20,100 ; f = -20,0,+20 %)
K_anchor = [1 5 20 100];
f_anchor = [-20 0 20];
fprintf('\nTable 4 anchor values (%% overshoot):\n');
fprintf('%8s %10s %10s %10s\n', 'K', 'f=-20%', 'f=0%', 'f=+20%');
for K = K_anchor
    fprintf('%8g %10.2f %10.2f %10.2f\n', K, PO(K,-20), PO(K,0), PO(K,20));
end

%% 7.7 Settling-time invariance and rise/peak-time trends
K_curve = logspace(log10(0.15), log10(300), 300);
zeta_K = zeta_cl(K_curve, 0);
wn_K   = wc*sqrt(1+K_curve);
wd_K   = wn_K.*sqrt(1-zeta_K.^2);
pole_real = -zeta_K.*wn_K;             % real part of the closed-loop poles -- should be constant

tp_K  = pi ./ wd_K;                                          % peak time
tr_K  = (pi - acos(zeta_K)) ./ wd_K;                         % 0-100% rise time
ts2_K = (log(50) - 0.5*log(1-zeta_K.^2)) ./ (zeta_K.*wn_K);  % 2% settling time (precise)

fprintf('\nClosed-loop pole real part: min=%.6f max=%.6f (target -wc/sqrt(2)=%.6f)\n', ...
        min(pole_real), max(pole_real), -wc/sqrt(2));
fprintf('Settling time (2%%): min=%.3fs max=%.3fs\n', min(ts2_K), max(ts2_K));
fprintf('Peak time: min=%.3fs max=%.3fs\n', min(tp_K), max(tp_K));

figure('Name','Timing Metrics and Pole Invariance');
subplot(1,2,1); hold on;
semilogx(K_curve, ts2_K, 'DisplayName', 'settling time (2%)');
semilogx(K_curve, tp_K, 'DisplayName', 'peak time');
semilogx(K_curve, tr_K, 'DisplayName', 'rise time (0-100%)');
grid on; legend show; xlabel('K (log scale)'); ylabel('Time (s)');
title('Timing Metrics vs. K');

subplot(1,2,2); hold on;
for K = [0.15 1 5 20 100 300]
    z = zeta_cl(K,0); wn = wc*sqrt(1+K); wd = wn*sqrt(1-z^2);
    plot(-z*wn, wd, 'o', 'MarkerSize', 8, 'DisplayName', sprintf('K=%g', K));
    plot(-z*wn, -wd, 'o', 'MarkerSize', 8, 'HandleVisibility', 'off');
end
xline(-wc/sqrt(2), '--', 'DisplayName', 'Re(s) = -\omega_c/\surd2');
grid on; legend show; xlabel('Real part'); ylabel('Imaginary part');
title('Closed-Loop Pole Locations vs. K (vertical locus)');

%% 7.8 Sensitivity structure and worst-case robust design
dPOdf = @(K,f) -pi*PO(K,f)./(1 - zeta_cl(K,f).^2).^1.5 .* (zeta0./(100*sqrt(1+K)));
sens_K = dPOdf(K_curve, 0);              % sensitivity evaluated at the nominal f = 0
[worst_sens, idx] = min(sens_K);         % most negative = largest sensitivity magnitude
K_worst_sens = K_curve(idx);
fprintf('\nMost sensitive point: K=%.3f, dPO/df=%.4f %%/point, PO there=%.2f%%\n', ...
        K_worst_sens, worst_sens, PO(K_worst_sens,0));

Kmax = @(PO_t,f) (zeta0*(1+f/100)).^2 .* (1 + pi^2./(log(100./PO_t)).^2) - 1;
PO_targets = [20 30 50 70];
K_nominal  = Kmax(PO_targets, 0);        % allowed gain, ignoring filter tolerance
K_robust   = Kmax(PO_targets, -30);      % allowed gain, guaranteeing the spec at f=-30%

fprintf('\n%8s %14s %14s %14s\n', 'PO target', 'K_max nominal', 'K_max robust', 'overdesign %%');
for i = 1:length(PO_targets)
    risk = (K_nominal(i)-K_robust(i))/K_nominal(i)*100;
    fprintf('%8g%% %14.2f %14.2f %14.1f\n', PO_targets(i), K_nominal(i), K_robust(i), risk);
end

figure('Name','Sensitivity and Robust Design');
subplot(1,2,1);
semilogx(K_curve, sens_K, 'LineWidth', 1.5); hold on;
plot(K_worst_sens, worst_sens, 'ko', 'MarkerSize', 8);
grid on; xlabel('K (log scale)'); ylabel('\partial PO/\partial f at f=0');
title('Overshoot Sensitivity vs. K');

subplot(1,2,2); hold on;
semilogx(K_curve, PO(K_curve,-30), 'DisplayName', 'worst case, f=-30%');
semilogx(K_curve, PO(K_curve,0),   'DisplayName', 'nominal, f=0%');
semilogx(K_curve, PO(K_curve,30),  'DisplayName', 'best case, f=+30%');
grid on; legend show; xlabel('K (log scale)'); ylabel('Percent overshoot (%)');
title('Worst/Nominal/Best-Case Overshoot Bounds');
