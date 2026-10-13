%% Discrete-Time Root Locus — Ts = 1 s
clear; clc; close all;

Ts = 1;

% Given discrete-time plant
z = tf('z', Ts);
Gz = 0.32*(z + 0.62)/((z - 1)*(z - 0.15));

% Create figure
figure('Color', 'w', ...
    'Name', 'Discrete Root Locus Ts = 1 s', ...
    'Position', [100 100 750 650]);

% Calculate root locus data
[r, k] = rlocus(Gz);

hold on;

% Plot each root-locus branch separately
branchColors = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980];

for i = 1:size(r, 1)
    plot(real(r(i,:)), imag(r(i,:)), ...
        'Color', branchColors(i,:), ...
        'LineWidth', 2.5, ...
        'DisplayName', sprintf('Root locus branch %d', i));
end

% Mark open-loop poles
poles = pole(Gz);
plot(real(poles), imag(poles), 'x', ...
    'Color', [0.85 0.1 0.1], ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'DisplayName', 'Open-loop poles');

% Mark open-loop zeros
zerosGz = zero(Gz);
plot(real(zerosGz), imag(zerosGz), 'o', ...
    'Color', [0.1 0.3 0.8], ...
    'MarkerSize', 10, ...
    'LineWidth', 2, ...
    'DisplayName', 'Open-loop zeros');

% Unit circle: discrete-time stability boundary
theta = linspace(0, 2*pi, 600);
plot(cos(theta), sin(theta), 'k--', ...
    'LineWidth', 2, ...
    'DisplayName', 'Unit circle');

% Formatting
axis equal;
xlim([-3.5 1.5]);
ylim([-1.5 1.5]);
grid on;
box on;

xlabel('Real(z)');
ylabel('Imaginary(z)');
title('Discrete-Time Root Locus (T_s = 1 s)');

ax = gca;
set(ax, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 15, ...
    'LineWidth', 3);

legend('Location', 'best', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 10);

% Export for report
if ~exist('assets', 'dir')
    mkdir('assets');
end

exportgraphics(gcf, ...
    'assets/root_locus_Ts_1s.png', ...
    'Resolution', 300);