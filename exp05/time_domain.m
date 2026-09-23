sys_def
Kp_vec = linspace(0.1, 30, 100);
Kp_valid=[];
for i = 1:length(Kp_vec)
    Kp = Kp_vec(i);
    
    % Closed-loop transfer function with unity feedback
    T = feedback(Kp * G, 1);
    
    % 1. Check Steady-State Error for a unit step (Target: <= 5% or 0.05)
    y_final = dcgain(T);
    err_ss = abs(1 - y_final);
    
    % 2 & 3. Check Transient Response (5% settling time <= 10s, Overshoot <= 2%)
    info = stepinfo(T, 'SettlingTimeThreshold', 0.05);
    settling_time = info.SettlingTime;
    overshoot = info.Overshoot;
    
    % Display results that meet all criteria
    if err_ss <= 0.05 && settling_time <= 10 && overshoot <= 2
        fprintf('Valid Kp found: %.2f (Error: %.2f%%, Settling Time: %.2fs, Overshoot: %.2f%%)\n', ...
            Kp, err_ss*100, settling_time, overshoot);
        Kp_valid(end+1)=Kp;
    end
end

Kp = 12.18;
T = feedback(Kp * G, 1);
[Gm,Pm] = margin(G);
GmdB = 20*log10(Gm)   % gain margin in dB
Pm  % phase margin in degrees
