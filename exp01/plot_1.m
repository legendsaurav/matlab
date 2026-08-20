clc;
clear;
close all;

%% =========================================================
% Settings
% =========================================================

K_values = [0.1 1 2 2.5 5 10 50];

% Create plots folder next to this MATLAB file
script_dir = fileparts(mfilename('fullpath'));

if isempty(script_dir)
    script_dir = pwd;
end

plot_dir = fullfile(script_dir, 'plots');

if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

%% =========================================================
% Generate and save plots
% =========================================================

for K = K_values

    % -----------------------------------------------------
    % Actual fourth-order system
    % -----------------------------------------------------

    G = zpk([-1 -0.2], ...
            [0 0 -1+3i -1-3i], ...
            5*K);

    T = feedback(G, 1);

    % -----------------------------------------------------
    % Second-order approximation
    % -----------------------------------------------------

    T2 = getT2(K);

    % -----------------------------------------------------
    % Create figure
    % -----------------------------------------------------

    fig = figure( ...
        'Color', 'white', ...
        'Position', [100 100 900 600]);

    % Plot actual system
    step(T);
    hold on;

    % Plot second-order approximation
    step(T2);

    % -----------------------------------------------------
    % Formatting
    % -----------------------------------------------------

    grid on;
    box on;

    ax = gca;
    ax.Color = 'white';
    ax.FontSize = 12;
    ax.LineWidth = 1;

    xlabel('Time (s)', ...
        'FontSize', 13);

    ylabel('Amplitude', ...
        'FontSize', 13);

    title(sprintf( ...
        'Step Response Comparison, K = %.2g', K), ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    legend( ...
        'Fourth-order system', ...
        'Second-order approximation', ...
        'Location', 'best', ...
        'FontSize', 11);

    % -----------------------------------------------------
    % Save figure
    % -----------------------------------------------------

    filename = fullfile( ...
        plot_dir, ...
        sprintf('step_response_K_%g.png', K));

    exportgraphics(fig, filename, ...
        'Resolution', 300, ...
        'BackgroundColor', 'white');

    close(fig);

    fprintf('Saved: %s\n', filename);

end

fprintf('\nAll plots saved to:\n%s\n', plot_dir);


%% =========================================================
% Second-order approximation function
% =========================================================

function T2 = getT2(K)

    % Open-loop transfer function
    G = zpk([-1 -0.2], ...
            [0 0 -1+3i -1-3i], ...
            5*K);

    % Closed-loop transfer function
    T = feedback(G, 1);

    % Closed-loop poles
    p = pole(T);

    % Select complex-conjugate poles
    p = p(imag(p) > 1e-8);

    % Select pair closest to imaginary axis
    [~, i] = max(real(p));
    pd = p(i);

    % Pole parameters
    sigma = -real(pd);
    wd = abs(imag(pd));

    wn = sqrt(sigma^2 + wd^2);
    zeta = sigma / wn;

    % Second-order approximation
    T2 = tf(wn^2, ...
            [1 2*zeta*wn wn^2]);

end