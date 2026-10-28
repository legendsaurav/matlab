%% EE208 Experiment 9 -- Core Analysis (Sections 5-6 of the report)
% Controller Design by Discrete Frequency Response
% Group 6B. Reproduces: Table 1-3, Figures 2-5.
% Requires base MATLAB + Control System Toolbox.

clear; clc; close all;

%% 5.1 Transfer function and pole-zero analysis
num_GOL = [1 1.7 -1];      % numerator of G_OL(z):  z^2 + 1.7z - 1
den_GOL = [1 0 0 -1];      % denominator of G_OL(z): z^3 - 1

zz = roots(num_GOL);       % zeros of G_OL(z)
pp = roots(den_GOL);       % poles of G_OL(z)

disp('Zeros of G_OL(z):'); disp(zz);
disp('Poles of G_OL(z):'); disp(pp);

figure('Name','Pole-Zero Map');
th = linspace(0,2*pi,400);
plot(cos(th), sin(th), 'k--'); hold on; grid on; axis equal;
plot(real(zz), imag(zz), 'o', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(pp), imag(pp), 'x', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Real Part'); ylabel('Imaginary Part'); title('Pole-Zero Map of G_{OL}(z)');
legend('Unit Circle','Zeros','Poles');

%% 5.2 Main-band frequency response
theta   = linspace(1e-3, pi-1e-4, 4000);
z_theta = exp(1j*theta);

Gz = polyval(num_GOL, z_theta) ./ polyval(den_GOL, z_theta);

magdB = 20*log10(abs(Gz));
phase = unwrap(angle(Gz))*180/pi;

figure('Name','Main-Band Response');
subplot(2,1,1); plot(theta, magdB); grid on; ylabel('Magnitude (dB)');
title('Main-Band Frequency Response of G_{OL}(z)');
subplot(2,1,2); plot(theta, phase); grid on; ylabel('Phase (deg)');
xlabel('Digital frequency \theta = \omega T_s (rad)');

%% 5.3 Boundary value and sampling-time invariance
G_at_boundary = polyval(num_GOL,-1) / polyval(den_GOL,-1);  % G_OL(e^{j*pi})
fprintf('G_OL(-1) = %.4f  ->  %.2f dB\n', G_at_boundary, 20*log10(abs(G_at_boundary)));

Ts_list = [1.0 0.5 0.2 0.1];   % sampling times under test
omegaN  = pi ./ Ts_list;        % corresponding Nyquist frequencies, rad/s

figure('Name','Effect of Ts');
for i = 1:length(Ts_list)
    Ts = Ts_list(i);
    omega = theta/Ts;
    subplot(2,2,i);
    plot(omega, magdB); grid on;
    xline(pi/Ts, '--');
    title(sprintf('T_s = %.1f s  (\\omega_N = %.2f rad/s)', Ts, pi/Ts));
    xlabel('\omega (rad/s)'); ylabel('Magnitude (dB)');
end

%% 5.4 Anti-aliasing filter design (Tustin, pre-warped)
theta_c   = pi/3;
Ts_design = 1.0;

wc = (2/Ts_design)*tan(theta_c/2);      % pre-warped analogue corner freq.
fprintf('wc = %.4f rad/s\n', wc);

C1_s = tf(wc, [1 wc]);                  % first-order:  C1(s) = wc /(s+wc)
C2_s = tf(wc^2, [1 sqrt(2)*wc wc^2]);   % 2nd-order Butterworth

C1_z = c2d(C1_s, Ts_design, 'tustin');
C2_z = c2d(C2_s, Ts_design, 'tustin');

[numC1, denC1] = tfdata(C1_z, 'v');
[numC2, denC2] = tfdata(C2_z, 'v');

disp('C2(z) numerator:');   disp(numC2);
disp('C2(z) denominator:'); disp(denC2);

H1 = polyval(numC1, z_theta) ./ polyval(denC1, z_theta);
H2 = polyval(numC2, z_theta) ./ polyval(denC2, z_theta);

casc1_dB = 20*log10(abs(H1.*Gz));
casc2_dB = 20*log10(abs(H2.*Gz));

figure('Name','Anti-Aliasing Filter');
subplot(2,1,1);
plot(theta, magdB, theta, 20*log10(abs(H1)), '--', theta, 20*log10(abs(H2)), '--');
grid on; legend('G_{OL}','C_1(z)','C_2(z)'); ylabel('Magnitude (dB)');
subplot(2,1,2);
plot(theta, magdB, theta, casc1_dB, theta, casc2_dB);
grid on; legend('Uncompensated','Compensated (order 1)','Compensated (order 2)');
xlabel('Digital frequency \theta = \omega T_s (rad)'); ylabel('Magnitude (dB)');

%% 5.5 Boundary and sideband evaluation
z_edge = exp(1j*(pi-1e-4));

G_edge  = polyval(num_GOL,z_edge) / polyval(den_GOL,z_edge);
H1_edge = polyval(numC1,z_edge)   / polyval(denC1,z_edge);
H2_edge = polyval(numC2,z_edge)   / polyval(denC2,z_edge);

fprintf('Boundary magnitude, uncompensated: %.2f dB\n', 20*log10(abs(G_edge)));
fprintf('Boundary magnitude, 1st-order AA : %.2f dB\n', 20*log10(abs(H1_edge*G_edge)));
fprintf('Boundary magnitude, 2nd-order AA : %.2f dB\n', 20*log10(abs(H2_edge*G_edge)));

omega_wide = linspace(-3.2*pi, 3.2*pi, 14000);   % Ts = 1 s, so omega = theta
z_wide = exp(1j*omega_wide);
Gw  = polyval(num_GOL, z_wide) ./ polyval(den_GOL, z_wide);
H2w = polyval(numC2, z_wide) ./ polyval(denC2, z_wide);

figure('Name','Sideband Response');
subplot(2,1,1); plot(omega_wide, 20*log10(abs(Gw))); grid on;
title('Uncompensated (sideband structure)'); ylabel('Magnitude (dB)');
subplot(2,1,2); plot(omega_wide, 20*log10(abs(H2w.*Gw))); grid on;
title('Compensated with 2nd-order AA filter'); ylabel('Magnitude (dB)');
xlabel('Physical frequency \omega (rad/s)');
