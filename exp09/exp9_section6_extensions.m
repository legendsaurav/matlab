%% EE208 Experiment 9 -- Section 6.4-6.7 Extensions
% (6.4) Impact of proportional gain k on the discrete frequency response
% (6.5) Nyquist plot of the compensated and uncompensated loop
% (6.6) Corner-frequency trade-off study
% (6.7) Time-domain aliasing demonstration
%
% Run exp9_core_analysis.m first (or paste its body above this file) --
% this script assumes num_GOL, den_GOL, theta, z_theta, Gz, magdB, Ts_design,
% wc, numC2, denC2, z_edge, G_edge all already exist, exactly as produced there.
% Requires base MATLAB + Control System Toolbox.

%% 6.4 Impact of proportional gain k
k_list = [0.1 1 5 10];

resPeak_dB  = zeros(size(k_list));   % resonance-peak magnitude for each k
edgePhase   = zeros(size(k_list));   % band-edge phase for each k
boundary_dB = zeros(size(k_list));   % boundary magnitude for each k

for i = 1:length(k_list)
    k = k_list(i);
    magdB_k = 20*log10(abs(k*Gz));
    resPeak_dB(i)  = max(magdB_k);
    Ledge          = k*polyval(num_GOL,z_edge)/polyval(den_GOL,z_edge);
    boundary_dB(i) = 20*log10(abs(Ledge));
    edgePhase(i)   = angle(Ledge)*180/pi;
end

figure('Name','Gain Sweep - Magnitude');
hold on;
for i = 1:length(k_list)
    plot(theta, 20*log10(abs(k_list(i)*Gz)), 'DisplayName', sprintf('k = %.1f', k_list(i)));
end
grid on; legend show; xlabel('\theta (rad)'); ylabel('Magnitude (dB)');
title('Effect of Proportional Gain k on 20log_{10}|k G_{OL}(e^{j\theta})|');

figure('Name','Gain Sweep - Peak and Phase vs k');
subplot(1,2,1); semilogx(k_list, resPeak_dB, 'o-'); grid on;
xlabel('k (log scale)'); ylabel('Resonance-peak magnitude (dB)');
subplot(1,2,2); semilogx(k_list, edgePhase, 'o-'); grid on; ylim([-5 5]);
xlabel('k (log scale)'); ylabel('Band-edge phase (deg)');

fprintf('Resonance peak (dB):  '); disp(resPeak_dB);
fprintf('Boundary value (dB):  '); disp(boundary_dB);
fprintf('Band-edge phase(deg): '); disp(edgePhase);   % should be ~identical for every k

%% 6.5 Nyquist plot, compensated vs uncompensated
theta_nyq = linspace(-pi, pi, 9001) + 1e-5;   % tiny offset avoids the 3 exact poles
z_nyq = exp(1j*theta_nyq);

Gz_nyq = polyval(num_GOL, z_nyq) ./ polyval(den_GOL, z_nyq);
C2_nyq = polyval(numC2,   z_nyq) ./ polyval(denC2,   z_nyq);

k_nyq_list = [0.5 2 5];

figure('Name','Nyquist Diagram');
subplot(1,2,1); hold on;
for k = k_nyq_list
    L = k*Gz_nyq;
    plot(real(L), imag(L), 'DisplayName', sprintf('k = %.1f', k));
end
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', '-1 point');
grid on; axis equal; xlim([-8 8]); ylim([-8 8]); legend show;
title('Uncompensated L(z) = k \cdot G_{OL}(z)'); xlabel('Real'); ylabel('Imag');

subplot(1,2,2); hold on;
for k = k_nyq_list
    L = k*C2_nyq.*Gz_nyq;
    plot(real(L), imag(L), 'DisplayName', sprintf('k = %.1f', k));
end
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', '-1 point');
grid on; axis equal; xlim([-8 8]); ylim([-8 8]); legend show;
title('Compensated L(z) = k \cdot C_2(z) \cdot G_{OL}(z)'); xlabel('Real'); ylabel('Imag');

%% 6.6 Corner-frequency trade-off study
theta_c_list = [pi/6, pi/3, pi/2, 2*pi/3];
boundary_dB_tc = zeros(size(theta_c_list));
resAtten_dB    = zeros(size(theta_c_list));
[~, idx_res]   = min(abs(theta - 2*pi/3));

figure('Name','Corner-Frequency Trade-off'); hold on;
plot(theta, magdB, 'k--', 'DisplayName', 'G_{OL} (uncompensated)');
for i = 1:length(theta_c_list)
    wc_i = (2/Ts_design)*tan(theta_c_list(i)/2);
    Ci_s = tf(wc_i^2, [1 sqrt(2)*wc_i wc_i^2]);
    Ci_z = c2d(Ci_s, Ts_design, 'tustin');
    [numCi, denCi] = tfdata(Ci_z, 'v');

    Hi = polyval(numCi, z_theta) ./ polyval(denCi, z_theta);
    casc_dB = 20*log10(abs(Hi.*Gz));
    plot(theta, casc_dB, 'DisplayName', sprintf('\\theta_c = %.4f rad', theta_c_list(i)));

    Hi_edge = polyval(numCi,z_edge) / polyval(denCi,z_edge);
    boundary_dB_tc(i) = 20*log10(abs(Hi_edge*G_edge));
    resAtten_dB(i)    = casc_dB(idx_res) - magdB(idx_res);
end
grid on; legend show; xlabel('\theta (rad)'); ylabel('Magnitude (dB)');
title('Cascaded Response for Different Corner Frequencies \theta_c');

fprintf('\ntheta_c (rad)   boundary (dB)   resonance attenuation (dB)\n');
for i = 1:length(theta_c_list)
    fprintf('%12.4f   %12.2f   %12.2f\n', theta_c_list(i), boundary_dB_tc(i), resAtten_dB(i));
end

%% 6.7 Time-domain aliasing demonstration
Ts = 1.0;
omega       = 1.2*pi/Ts;         % physical freq, just above Nyquist pi/Ts
omega_alias = 2*pi/Ts - omega;   % fold-back frequency

n = 0:20;
x_n = sin(omega*n*Ts);                    % what the ADC actually records
x_alias_check = -sin(omega_alias*n*Ts);   % sign-corrected fold-back (see report note)
max_err = max(abs(x_n - x_alias_check));  % ~0 confirms the fold-back identity
fprintf('\nmax|x_n - alias_check| = %.3e  (confirms exact fold-back identity)\n', max_err);

s_val  = 1j*omega;
C2_val = wc^2 / (s_val^2 + sqrt(2)*wc*s_val + wc^2);   % steady-state filter response at omega
mag_C2 = abs(C2_val);  ph_C2 = angle(C2_val);
y_n = mag_C2*sin(omega*n*Ts + ph_C2);      % samples after the analogue pre-filter
fprintf('|C2(j*omega)| = %.4f (%.2f dB), phase = %.2f deg\n', mag_C2, 20*log10(mag_C2), ph_C2*180/pi);

t_fine = linspace(0, 20, 4000);
figure('Name','Time-Domain Aliasing');
subplot(2,1,1); hold on;
plot(t_fine, sin(omega*t_fine), 'Color', [0.7 0.7 0.7]);
plot(t_fine, -sin(omega_alias*t_fine), 'r', 'LineWidth', 1.5);
stem(n*Ts, x_n, 'filled');
grid on; ylim([-1.3 1.3]); title('Direct Sampling (No Pre-Filter): Aliasing Occurs');
legend('true sin(\omega t)', 'alias -sin(\omega_{alias}t)', 'samples x[n]');

subplot(2,1,2); hold on;
plot(t_fine, mag_C2*sin(omega*t_fine + ph_C2), 'g', 'LineWidth', 1.5);
stem(n*Ts, y_n, 'filled');
grid on; ylim([-1.3 1.3]); title('With Analog Pre-Filter Before Sampling: Alias Suppressed');
legend('filtered C_2(s) output', 'samples y[n]'); xlabel('Time (s)');
