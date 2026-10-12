%% VERIFY_REPORT_NUMBERS  Pass/fail check of the function library.
% Plant (lab sheet, T = 1 s):  G(z) = 0.32 (z + 0.62) / [ (z-1)(z-0.15) ]
% Run this first. It needs only the nine function files (no toolbox except
% for the c2d check, which uses the Control System Toolbox).
%
% Checks
%   1. Kp,min / Kp,max / Kp,stab for T = 0.1 ... 1.0 s on the recovered
%      continuous model G(s) = K/[s(s+a)] against reference values that were
%      computed independently (Python/SciPy).
%   2. The same three gains for the EXACT lab-sheet coefficients at T = 1 s.
%   3. zoh_coeffs against MATLAB's own c2d(...,'zoh').
clear; clc;

pz_given = 0.15;  b1_given = 0.32;  zero_given = 0.62;  T1 = 1;
b0_given = b1_given*zero_given;

a = -log(pz_given)/T1;                              % 1.8971
K = b1_given*a^2/(a*T1 - 1 + pz_given);             % 1.0999
zeta_min = -log(0.05)/sqrt(pi^2 + log(0.05)^2);     % 0.6901

Ts      = (1:10)/10;
ref_min = [0.7803 0.7447 0.7111 0.6794 0.6497 0.6217 0.5955 0.5709 0.5478 0.5261];
ref_max = [1.5613 1.4299 1.3182 1.2224 1.1395 1.0671 1.0035 0.9473 0.8972 0.8524];
ref_stb = [35.623 18.410 12.697 9.859 8.172 7.061 6.279 5.702 5.262 4.919];

tol = 5e-4;  nfail = 0;
fprintf('--- 1. T sweep on the recovered continuous model ---\n');
fprintf('  T     Kp,min   (ref)    Kp,max   (ref)    Kp,stab  (ref)\n');
for i = 1:10
    [b1,b0,p] = zoh_coeffs(Ts(i),a,K);
    [kmin,kbi] = breakaway_gain(b1,b0,p);
    kmax = zeta_crossing_gain(b1,b0,p,zeta_min,kmin*1.0005,kbi*0.999);
    kstb = marginal_stability_gain(b0,p);
    ok = abs(kmin-ref_min(i))<tol && abs(kmax-ref_max(i))<tol && abs(kstb-ref_stb(i))<2e-3;
    nfail = nfail + ~ok;
    fprintf('%4.1f  %7.4f  %7.4f  %7.4f  %7.4f  %7.3f  %7.3f   %s\n', ...
        Ts(i),kmin,ref_min(i),kmax,ref_max(i),kstb,ref_stb(i),ternary(ok,'PASS','FAIL'));
end

fprintf('\n--- 2. Exact lab-sheet G(z) at T = 1 s ---\n');
[kmin,kbi] = breakaway_gain(b1_given,b0_given,pz_given);
kmax = zeta_crossing_gain(b1_given,b0_given,pz_given,zeta_min,kmin*1.0005,kbi*0.999);
kstb = marginal_stability_gain(b0_given,pz_given);
ref = [0.4883 0.7845 4.284];
ok2 = abs(kmin-ref(1))<tol && abs(kmax-ref(2))<tol && abs(kstb-ref(3))<2e-3;
nfail = nfail + ~ok2;
fprintf('Kp,min = %.4f (ref %.4f)   Kp,max = %.4f (ref %.4f)   Kp,stab = %.3f (ref %.3f)   %s\n', ...
    kmin,ref(1),kmax,ref(2),kstb,ref(3),ternary(ok2,'PASS','FAIL'));

fprintf('\n--- 3. Closed-form ZOH vs c2d ---\n');
Gz  = c2d(tf(K,[1 a 0]),1,'zoh');
[b1,b0,p] = zoh_coeffs(1,a,K);
okc2d = max(abs(Gz.Numerator{1}-[0 b1 b0])) < 1e-9 && abs(Gz.Denominator{1}(3)-p) < 1e-9;
nfail = nfail + ~okc2d;
fprintf('zoh_coeffs vs c2d(...,''zoh'') at T=1: %s\n', ternary(okc2d,'PASS','FAIL'));

fprintf('\nNote: the zero implied by the recovered ZOH model at T=1 s is -%.3f, the lab sheet gives -%.2f.\n', ...
    b0/b1, zero_given);

if nfail == 0, fprintf('\nALL CHECKS PASSED\n'); else, fprintf('\n%d CHECK(S) FAILED\n',nfail); end

function out = ternary(c,x,y)
    if c, out = x; else, out = y; end
end
