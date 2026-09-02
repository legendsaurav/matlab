%% ========================================================================
%  MAIN ENTRY POINT FOR EXPERIMENT 3
%  Runs the full MATLAB analysis workflow in the correct order and writes
%  generated figures into the project asset folders.
%  This script also keeps the output analysis organized under:
%    - assets/ for generated figures
%    - artifacts/ for exported numeric tables and CSV summaries
%    - report/ for the final report files
% ========================================================================
clc; clear; close all;

baseDir = fileparts(mfilename('fullpath'));
assetsDir = fullfile(baseDir, 'assets');
artifactsDir = fullfile(baseDir, 'artifacts');
reportDir = fullfile(baseDir, 'report');

if ~exist(assetsDir, 'dir'), mkdir(assetsDir); end
if ~exist(artifactsDir, 'dir'), mkdir(artifactsDir); end
if ~exist(reportDir, 'dir'), mkdir(reportDir); end

fprintf('=== Experiment 3 main workflow started ===\n');
fprintf('Assets folder: %s\n', assetsDir);
fprintf('Artifacts folder: %s\n', artifactsDir);
fprintf('Report folder: %s\n', reportDir);

run(fullfile(baseDir, 'matlab_partA_rootlocus.m'));
run(fullfile(baseDir, 'matlab_partB_damping_stability_sensitivity.m'));
run(fullfile(baseDir, 'matlab_partC_overshoot_settling_step.m'));
run(fullfile(baseDir, 'matlab_partD_extended_analysis.m'));
run(fullfile(baseDir, 'matlab_partE_multi_a_and_heatmaps.m'));

fprintf('\n=== Experiment 3 complete ===\n');
fprintf('Figures organized under %s\n', assetsDir);
fprintf('Numeric outputs are in %s\n', artifactsDir);
fprintf('Report documents remain in %s\n', reportDir);
