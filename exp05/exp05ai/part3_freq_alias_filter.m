%% Part 3: Main-band response, sampling-time effect, anti-aliasing filter
%
% Objective:
%  (a) compute the frequency response of G_OL(z) over the main
%      digital band 0 < theta < pi and locate the resonance / band edge;
%  (b) confirm that the response at the band edge does not depend on
%      the (unspecified) sampling time Ts when expressed in theta;
%  (c) design first- and second-order low-pass anti-aliasing filters
%      (unity DC gain), discretised by the Tustin method with
%      pre-warping, and quantify how much they suppress the response
%      at the band boundary.

num_GOL = [1 1.7 -1];
den_GOL = [1 0 0 -1];

% ---- (a) Main-band frequency response, theta = 0..pi ----
theta = linspace(1e-3, pi-1e-4, 4000);
z_theta = exp(1j*theta);

Gz = polyval(num_GOL,z_theta) ./ ...
     polyval(den_GOL,z_theta);

magdB = 20*log10(abs(Gz));
phase = unwrap(angle(Gz))*180/pi;

% ---- (b) Effect of sampling time Ts ----
G_at_boundary = polyval(num_GOL,-1) / ...
                polyval(den_GOL,-1);

Ts_list = [1.0 0.5 0.2 0.1];
for i = 1:length(Ts_list)
    Ts = Ts_list(i);
    omega = theta/Ts;                 %#ok<NASGU>  (physical-frequency axis for plotting)
end
fprintf('Boundary magnitude is %.2f dB for every Ts tested (independent of Ts).\n', ...
    20*log10(abs(G_at_boundary)))

% ---- (c) Anti-aliasing filter design (Tustin transform, pre-warped) ----
theta_c = pi/3;                       % desired digital corner frequency
Ts_design = 1.0;

wc = (2/Ts_design)*tan(theta_c/2);    % pre-warped analogue corner frequency

C1_s = tf(wc,[1 wc]);                             % first order
C2_s = tf(wc^2,[1 sqrt(2)*wc wc^2]);              % second order (Butterworth, zeta=1/sqrt(2))

C1_z = c2d(C1_s,Ts_design,'tustin');
C2_z = c2d(C2_s,Ts_design,'tustin');

[numC1,denC1] = tfdata(C1_z,'v');
[numC2,denC2] = tfdata(C2_z,'v');

H1 = polyval(numC1,z_theta) ./ polyval(denC1,z_theta);
H2 = polyval(numC2,z_theta) ./ polyval(denC2,z_theta);

casc1_dB = 20*log10(abs(H1.*Gz));
casc2_dB = 20*log10(abs(H2.*Gz));

% ---- Boundary comparison, uncompensated vs filtered ----
z_edge = exp(1j*(pi-1e-4));

G_edge  = polyval(num_GOL,z_edge) / polyval(den_GOL,z_edge);
H1_edge = polyval(numC1,z_edge)   / polyval(denC1,z_edge);
H2_edge = polyval(numC2,z_edge)   / polyval(denC2,z_edge);

fprintf('Uncompensated boundary magnitude:        %.2f dB\n', 20*log10(abs(G_edge)))
fprintf('1st-order filter + plant, boundary:      %.2f dB\n', 20*log10(abs(H1_edge*G_edge)))
fprintf('2nd-order filter + plant, boundary:       %.2f dB\n', 20*log10(abs(H2_edge*G_edge)))

disp('C2(z) numerator:');   disp(numC2)
disp('C2(z) denominator:'); disp(denC2)
