function [OS,ts] = evaluate_design(b1,b0,p,Kp,T)
%EVALUATE_DESIGN  True (simulated) percent overshoot and 2%-criterion
%settling time of the closed digital loop at one design point.
%
%   [OS,ts] = evaluate_design(b1,b0,p,Kp,T)
%
%   Builds the ZOH plant from (b1,b0,p), closes the loop with unity
%   feedback and gain Kp, and reads the full discrete step response. This
%   is the "simulated" (not formula-based) OS/ts used throughout the
%   report.
%
%   Requires: Control System Toolbox (tf, feedback, stepinfo).
%   See EE208 Experiment 7 report, Sec. 4.6.

    Gz   = tf([b1 b0],[1 -(1+p) p], T);
    Gcl  = feedback(Kp*Gz, 1);
    info = stepinfo(Gcl,'SettlingTimeThreshold',0.02);
    OS = info.Overshoot;
    ts = info.SettlingTime;
end
