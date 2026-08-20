clc;
clear;
close all;

%% =========================================================
% Settings
% =========================================================

K_values = [0.1 0.25 0.5 1 1.5 2 2.5 3 ...
            4 5 7.5 10 20 50 100];

% Folder containing this MATLAB file
script_dir = fileparts(mfilename('fullpath'));

if isempty(script_dir)
    script_dir = pwd;
end

% Output folders
artifact_dir = fullfile(script_dir, 'artifacts');
asset_dir    = fullfile(script_dir, 'assets');

if ~exist(artifact_dir, 'dir')
    mkdir(artifact_dir);
end

if ~exist(asset_dir, 'dir')
    mkdir(asset_dir);
end

%% =========================================================
% Results table
% =========================================================

results = table();

%% =========================================================
% Main loop
% =========================================================

for K = K_values

    %% -----------------------------------------------------
    % Nominal open-loop poles
    % -----------------------------------------------------

    p_nom_1 = -1 + 3i;
    p_nom_2 = -1 - 3i;

    %% -----------------------------------------------------
    % -20% pole variation
    % -----------------------------------------------------

    p_m20_1 = -0.8 + 2.4i;
    p_m20_2 = -0.8 - 2.4i;

    %% -----------------------------------------------------
    % +20% pole variation
    % -----------------------------------------------------

    p_p20_1 = -1.2 + 3.6i;
    p_p20_2 = -1.2 - 3.6i;

    %% -----------------------------------------------------
    % Open-loop transfer functions
    % -----------------------------------------------------

    G_nom = zpk([-1 -0.2], ...
                [0 0 p_nom_1 p_nom_2], ...
                5*K);

    G_m20 = zpk([-1 -0.2], ...
                [0 0 p_m20_1 p_m20_2], ...
                5*K);

    G_p20 = zpk([-1 -0.2], ...
                [0 0 p_p20_1 p_p20_2], ...
                5*K);

    %% -----------------------------------------------------
    % Closed-loop systems
    % -----------------------------------------------------

    T_nom = feedback(G_nom,1);
    T_m20 = feedback(G_m20,1);
    T_p20 = feedback(G_p20,1);

    %% -----------------------------------------------------
    % Step-response metrics
    % -----------------------------------------------------

    info_nom = stepinfo(T_nom);
    info_m20 = stepinfo(T_m20);
    info_p20 = stepinfo(T_p20);

    PO_nom  = info_nom.Overshoot;
    Tp_nom  = info_nom.PeakTime;
    Ts_nom  = info_nom.SettlingTime;
    ess_nom = abs(1 - dcgain(T_nom));

    PO_m20  = info_m20.Overshoot;
    Tp_m20  = info_m20.PeakTime;
    Ts_m20  = info_m20.SettlingTime;
    ess_m20 = abs(1 - dcgain(T_m20));

    PO_p20  = info_p20.Overshoot;
    Tp_p20  = info_p20.PeakTime;
    Ts_p20  = info_p20.SettlingTime;
    ess_p20 = abs(1 - dcgain(T_p20));

    %% -----------------------------------------------------
    % Percentage deviation from nominal
    % -----------------------------------------------------

    PO_error_m20 = abs(PO_m20 - PO_nom) / abs(PO_nom) * 100;
    Tp_error_m20 = abs(Tp_m20 - Tp_nom) / abs(Tp_nom) * 100;
    Ts_error_m20 = abs(Ts_m20 - Ts_nom) / abs(Ts_nom) * 100;

    PO_error_p20 = abs(PO_p20 - PO_nom) / abs(PO_nom) * 100;
    Tp_error_p20 = abs(Tp_p20 - Tp_nom) / abs(Tp_nom) * 100;
    Ts_error_p20 = abs(Ts_p20 - Ts_nom) / abs(Ts_nom) * 100;

    %% -----------------------------------------------------
    % Store results
    % -----------------------------------------------------

    newRow = table( ...
        K, ...
        PO_nom, Tp_nom, Ts_nom, ess_nom, ...
        PO_m20, Tp_m20, Ts_m20, ess_m20, ...
        PO_p20, Tp_p20, Ts_p20, ess_p20, ...
        PO_error_m20, Tp_error_m20, Ts_error_m20, ...
        PO_error_p20, Tp_error_p20, Ts_error_p20, ...
        'VariableNames', { ...
        'K', ...
        'Nominal_OS_percent', ...
        'Nominal_Tp', ...
        'Nominal_Ts', ...
        'Nominal_ess', ...
        'Minus20_OS_percent', ...
        'Minus20_Tp', ...
        'Minus20_Ts', ...
        'Minus20_ess', ...
        'Plus20_OS_percent', ...
        'Plus20_Tp', ...
        'Plus20_Ts', ...
        'Plus20_ess', ...
        'Minus20_OS_error_percent', ...
        'Minus20_Tp_error_percent', ...
        'Minus20_Ts_error_percent', ...
        'Plus20_OS_error_percent', ...
        'Plus20_Tp_error_percent', ...
        'Plus20_Ts_error_percent'});

    results = [results; newRow];

    %% =====================================================
    % Step-response plot
    % =====================================================

    fig = figure( ...
        'Color','white', ...
        'Position',[100 100 900 600]);

    step(T_nom);
    hold on;

    step(T_m20);
    step(T_p20);

    grid on;
    box on;

    ax = gca;

    xlabel('Time (s)', 'FontSize', 13);
    ylabel('Amplitude', 'FontSize', 13);

    title(sprintf( ...
        'Effect of \pm20%% Open-Loop Pole Variation, K = %.2g', K), ...
        'FontSize', 14, ...
        'FontWeight','bold');

    legend( ...
        'Nominal', ...
        '-20% pole variation', ...
        '+20% pole variation', ...
        'Location','best', ...
        'FontSize',11);

    filename = fullfile( ...
        asset_dir, ...
        sprintf('pole_variation_K_%g.png',K));

    exportgraphics(fig, filename, ...
        'Resolution',300, ...
        'BackgroundColor','white');

    close(fig);

    fprintf('K = %.2g completed\n',K);

end

%% =========================================================
% Write CSV
% =========================================================

csv_file = fullfile( ...
    artifact_dir, ...
    'pole_variation_results.csv');

writetable(results, csv_file);

fprintf('\n============================================\n');
fprintf('Branch B completed.\n');
fprintf('CSV saved to:\n%s\n',csv_file);
fprintf('Plots saved to:\n%s\n',asset_dir);
fprintf('============================================\n');