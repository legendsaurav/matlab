%% EE208 Control Systems Lab - Experiment 1

clear;
close all;
clc;

K_values = [0.1 0.25 0.5 1 1.5 2 2.5 3 4 5 7.5 10 20 50 100];

results = table();

for K = K_values

    % Open-loop and closed-loop transfer functions
    G = zpk([-1 -0.2], [0 0 -1+3i -1-3i], 5*K);
    T = feedback(G, 1);

    % -------------------------------------------------
    % Actual 4th-order system
    % -------------------------------------------------
    [PO,Tp,Ts,ess] = fourthOrderMetrics(T);

    actualPoles = pole(T);

    % -------------------------------------------------
    % Second-order approximation
    % -------------------------------------------------
    [PO2,Tp2,Ts2,ess2,secondOrderPoles,T2] = ...
        secondOrderMetrics(T);

    % -------------------------------------------------
    % Percentage errors
    % -------------------------------------------------
    OS_error = abs(PO2 - PO) / abs(PO) * 100;
    Tp_error = abs(Tp2 - Tp) / abs(Tp) * 100;
    Ts_error = abs(Ts2 - Ts) / abs(Ts) * 100;

    % -------------------------------------------------
    % Store results
    % -------------------------------------------------
    newRow = table( ...
        K, ...
        PO, Tp, Ts, ess, ...
        PO2, Tp2, Ts2, ess2, ...
        OS_error, Tp_error, Ts_error, ...
        actualPoles(1), actualPoles(2), ...
        actualPoles(3), actualPoles(4), ...
        secondOrderPoles(1), secondOrderPoles(2), ...
        'VariableNames', { ...
        'K', ...
        'Actual_OS_percent', ...
        'Actual_Tp', ...
        'Actual_Ts', ...
        'Actual_ess', ...
        'SecondOrder_OS_percent', ...
        'SecondOrder_Tp', ...
        'SecondOrder_Ts', ...
        'SecondOrder_ess', ...
        'OS_Error_percent', ...
        'Tp_Error_percent', ...
        'Ts_Error_percent', ...
        'Pole1', 'Pole2', 'Pole3', 'Pole4', ...
        'DominantPole1', 'DominantPole2'});

    results = [results; newRow];

end

% -------------------------------------------------
% Write results to CSV
% -------------------------------------------------

writetable(results, 'experiment1_results.csv');

disp('Results written to experiment1_results.csv');


%% =========================================================
% Actual 4th-order step-response metrics
% =========================================================
function [PO,Tp,Ts,ess] = fourthOrderMetrics(T)

info = stepinfo(T);

PO  = info.Overshoot;
Tp  = info.PeakTime;
Ts  = info.SettlingTime;
ess = abs(1 - dcgain(T));

end


%% =========================================================
% Second-order approximation metrics
% =========================================================
function [PO,Tp,Ts,ess,p2,T2] = secondOrderMetrics(T)

p = pole(T);

% Keep poles with positive imaginary part
p = p(imag(p) > 1e-8);

% Dominant complex pole = closest to imaginary axis
[~,i] = max(real(p));
pd = p(i);

% Dominant conjugate pair
p2 = [pd; conj(pd)];

% -------------------------------------------------
% Second-order parameters
% -------------------------------------------------

sigma = -real(pd);
wd = abs(imag(pd));

wn = sqrt(sigma^2 + wd^2);
zeta = sigma/wn;

% -------------------------------------------------
% Construct second-order approximation
%
% Unity DC gain is used because the original
% closed-loop system has unity steady-state gain.
% -------------------------------------------------

T2 = tf(wn^2, [1 2*zeta*wn wn^2]);

% -------------------------------------------------
% Step-response metrics
% -------------------------------------------------

info = stepinfo(T2);

PO  = info.Overshoot;
Tp  = info.PeakTime;
Ts  = info.SettlingTime;
ess = abs(1 - dcgain(T2));

end