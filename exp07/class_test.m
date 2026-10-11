clc;
clear;

% Given digital plant at Ts = 1 s
z = tf('z',1);
Gz = 0.32*(z+0.62)/((z-1)*(z-0.15));

% Recover continuous-time plant
Gs = d2c(Gz,'zoh');

Kp = 1;

% Continuous-time closed loop
Tc = feedback(Kp*Gs,1);
info_c = stepinfo(Tc);

fprintf('========== CONTINUOUS DOMAIN ==========\n');
fprintf('Overshoot      = %.2f %%\n',info_c.Overshoot);
fprintf('Settling Time  = %.2f s\n\n',info_c.SettlingTime);

% Different sampling times
T_grid = 0.1:0.1:1;

fprintf('========== DISCRETE DOMAIN ==========\n');

for sampleTime = T_grid

    % Correct discrete model for this sampling time
    Gz_new = c2d(Gs,sampleTime,'zoh');

    % Closed-loop system
    T = feedback(Kp*Gz_new,1);

    % Step information
    info = stepinfo(T);

    fprintf('Ts = %.1f s | OS = %.2f %% | Settling = %.2f s\n',...
        sampleTime,...
        info.Overshoot,...
        info.SettlingTime);
end
