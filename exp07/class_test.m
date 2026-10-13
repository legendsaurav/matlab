
%% Discretized Plants and Root Loci for Different Sampling Times
clear; clc; close all;

Kp = 1;
T_grid = 0.1:0.1:1.0;

%% 1. Reconstruct continuous-time plant from Ts = 1 s model
z = tf('z', 1);
Gz = 0.32*(z + 0.62)/((z - 1)*(z - 0.15));

Gs = d2c(Gz, 'zoh');

disp('Continuous-time equivalent plant:');
Gs

disp('Continuous-time plant zeros:');
disp(zero(Gs));

%% 2. Create table: sampling time versus discretized plant
Gd_list = cell(numel(T_grid), 1);
Gd_string = strings(numel(T_grid), 1);

for i = 1:numel(T_grid)
    Ts = T_grid(i);

    Gd_list{i} = c2d(Gs, Ts, 'zoh');

    [num, den] = tfdata(Gd_list{i}, 'v');
    Gd_string(i) = sprintf('(%s) / (%s)', ...
        strtrim(poly2str(num, 'z')), ...
        strtrim(poly2str(den, 'z')));
end

PlantTable = table(T_grid(:), Gd_string, ...
    'VariableNames', {'Ts_s', 'Gd_z'});

disp('TABLE 1: Sampling time versus discrete-time plant');
disp(PlantTable);

%% 3. Continuous-time reference response
Tc = feedback(Kp*Gs, 1);
info_c = stepinfo(Tc);

fprintf('\nCONTINUOUS-TIME REFERENCE (Kp = %.2f)\n', Kp);
fprintf('Overshoot     = %.2f %%\n', info_c.Overshoot);
fprintf('Settling time = %.3f s\n', info_c.SettlingTime);
disp('Closed-loop poles:');
disp(pole(Tc));

%% 4. Plot and export separate root-locus figures
if ~exist('assets', 'dir')
    mkdir('assets');
end

theta = linspace(0, 2*pi, 600);

branchColors = [
    0.0000 0.4470 0.7410
    0.8500 0.3250 0.0980
];

for i = 1:numel(T_grid)
    Ts = T_grid(i);
    Gd = Gd_list{i};

    % Calculate root-locus data
    [r, ~] = rlocus(Gd);

    % Create a separate figure
    fig = figure('Color', 'w', ...
        'Name', sprintf('Root Locus Ts = %.1f s', Ts), ...
        'Units', 'pixels', ...
        'Position', [100 100 1000 800]);

    hold on;

    % Plot each branch separately
    for j = 1:size(r, 1)
        plot(real(r(j,:)), imag(r(j,:)), ...
            'Color', branchColors(mod(j-1,2)+1,:), ...
            'LineWidth', 2.5, ...
            'DisplayName', sprintf('Root locus branch %d', j));
    end

    % Open-loop poles
    p = pole(Gd);
    plot(real(p), imag(p), 'rx', ...
        'MarkerSize', 12, ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Open-loop poles');

    % Open-loop zeros
    q = zero(Gd);
    plot(real(q), imag(q), 'bo', ...
        'MarkerSize', 10, ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Open-loop zeros');

    % Unit circle: stability boundary
    plot(cos(theta), sin(theta), 'k--', ...
        'LineWidth', 1.8, ...
        'DisplayName', 'Unit circle');

    % Formatting
    grid on;
    box on;

    xlabel('Real(z)');
    ylabel('Imaginary(z)');
    title(sprintf('Discrete-Time Root Locus (T_s = %.1f s)', Ts));

    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 15, ...
        'LineWidth', 1.5);

    legend('Location', 'best', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11);

    % Export each figure separately
    filename = sprintf('root_locus_Ts_%0.1fs', Ts);
    filename = strrep(filename, '.', 'p');

    exportgraphics(fig, ...
        fullfile('assets', [filename '.png']), ...
        'Resolution', 600, ...
        'BackgroundColor', 'white');

    exportgraphics(fig, ...
        fullfile('assets', [filename '.pdf']), ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

    close(fig);
end

fprintf('\nSaved separate root-locus PNG and PDF figures in assets/.\n');%% 4. Plot and export separate root-locus figures
if ~exist('assets', 'dir')
    mkdir('assets');
end

theta = linspace(0, 2*pi, 600);

branchColors = [
    0.0000 0.4470 0.7410
    0.8500 0.3250 0.0980
    ];

for i = 1:numel(T_grid)
    Ts = T_grid(i);
    Gd = Gd_list{i};

    % Calculate root-locus data
    [r, ~] = rlocus(Gd);

    % Wide figure with dedicated legend space
    fig = figure('Color', 'w', ...
        'Name', sprintf('Root Locus Ts = %.1f s', Ts), ...
        'Units', 'pixels', ...
        'Position', [100 100 1400 800]);

    ax = axes(fig, 'Position', [0.09 0.12 0.62 0.78]);
    hold(ax, 'on');

    % Plot each branch separately
    for j = 1:size(r, 1)
        plot(ax, real(r(j,:)), imag(r(j,:)), ...
            'Color', branchColors(mod(j-1,2)+1,:), ...
            'LineWidth', 2.5, ...
            'DisplayName', sprintf('Root locus branch %d', j));
    end

    % Open-loop poles
    p = pole(Gd);
    plot(ax, real(p), imag(p), 'rx', ...
        'MarkerSize', 12, ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Open-loop poles');

    % Open-loop zeros
    q = zero(Gd);
    plot(ax, real(q), imag(q), 'bo', ...
        'MarkerSize', 10, ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Open-loop zeros');

    % Unit circle
    plot(ax, cos(theta), sin(theta), 'k--', ...
        'LineWidth', 1.8, ...
        'DisplayName', 'Unit circle');

    % Scale limits relative to 1/Ts, retaining the unit circle
    scale = 1/Ts;
    xlim(ax, [-max(1.5, 3.5*scale) max(1.5, 1.5*scale)]);
    ylim(ax, [-max(1.5, 1.5*scale) max(1.5, 1.5*scale)]);

    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');

    xlabel(ax, 'Real(z)');
    ylabel(ax, 'Imaginary(z)');
    title(ax, sprintf('Discrete-Time Root Locus (T_s = %.1f s)', Ts));

    set(ax, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 15, ...
        'LineWidth', 1.5);

    % Dedicated legend area, aligned vertically
    lgd = legend(ax, 'show', ...
        'Location', 'none', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 12);

    lgd.Units = 'normalized';
    lgd.Position = [0.75 0.35 0.22 0.30];
    lgd.Box = 'on';

    % Export separately
    filename = strrep(sprintf('root_locus_Ts_%0.1fs', Ts), '.', 'p');

    exportgraphics(fig, ...
        fullfile('assets', [filename '.png']), ...
        'Resolution', 600, ...
        'BackgroundColor', 'white');

    exportgraphics(fig, ...
        fullfile('assets', [filename '.pdf']), ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

    close(fig);
end

fprintf('\nSaved separate root-locus PNG and PDF figures in assets/.\n');

%% 5. Export figure
if ~exist('assets', 'dir')
    mkdir('assets');
end

exportgraphics(fig, ...
    'assets/discrete_root_loci_all_Ts.png', ...
    'Resolution', 600);

exportgraphics(fig, ...
    'assets/discrete_root_loci_all_Ts.pdf', ...
    'ContentType', 'vector');