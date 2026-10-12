%% EE208 Experiment 7 -- Controller Design on MATLAB Using Discrete Root Loci
% MAIN SCRIPT.  Plant for this version:
%
%                    0.32 (z + 0.62)
%      G(z)|T=1s = -----------------------        (handwritten lab sheet)
%                   (z - 1)(z - 0.15)
%
% Reproduces the whole experiment end to end:
%   Part I  : base proportional design (Objectives 1-5)
%   Part II : extended analysis: graphical workflow, margins, Nyquist,
%             bandwidth, zero migration, sensitivity, robustness, lead
%             compensator, nomograph, intersample peak.
%
% TO USE A DIFFERENT PLANT edit ONLY the "0. SETTINGS" block below.
% No design gain is hard-coded anywhere; every Kp is computed.
%
% HOW TO RUN
%   1. Keep this file and the nine function files in ONE folder:
%        zoh_coeffs  effective_zeta  breakaway_gain  zeta_error
%        zeta_crossing_gain  evaluate_design  marginal_stability_gain
%        fit_inverse_law  pole_sensitivity
%   2. cd to that folder, run verify_report_numbers, then run this script.
%
% REQUIREMENTS
%   MATLAB R2018b+ (yline/xline; script-local functions need R2016b+).
%   Control System Toolbox: tf c2d ss feedback stepinfo margin bandwidth
%   nyquistplot rlocus rlocusplot zgrid.  fzero/fminsearch are base MATLAB.

clear; close all; clc;

%% 0. SETTINGS  (edit here for a different plant / spec)
T1         = 1;        % sampling time at which G(z) is given (s)
pz_given   = 0.15;     % open-loop finite pole  (z - 0.15)
b1_given   = 0.32;     % numerator gain         0.32 (z + 0.62)
zero_given = 0.62;     % numerator zero at z = -0.62
OS_spec    = 5;        % percent overshoot limit
ts_spec    = 9;        % 2% settling-time limit (s)
ts_lead    = 5;        % settling-time target for the lead compensator (s); must beat the best proportional-only ts (printed in Sec. 11)
lead_margin= 1;        % required margin under ts_lead (s); samples are 1 s apart
OS_tol     = 0.05;     % percentage-point tolerance when testing OS <= spec

% Which plant description do the T = 1 s analyses use?
%   true  : the exact G(z) written on the lab sheet (b0 = b1*zero_given)
%   false : the ZOH discretisation of the recovered continuous plant K/[s(s+a)]
% The two differ slightly (see the mismatch report printed in Sec. 4.1).
% The T-sweep (T = 0.1..1.0) always needs the continuous model.
USE_GIVEN_AT_T1 = true;

%% 4.1 Recover G(s) from G(z)  (Objective 1)
b0_given = b1_given*zero_given;
a  = -log(pz_given)/T1;                          % z = exp(-aT)
p1 = exp(-a*T1);
K  = b1_given*a^2/(a*T1 - 1 + p1);               % match the z^1 coefficient

Gs       = tf(K,[1 a 0]);                        % recovered plant K/[s(s+a)]
Gz_check = c2d(Gs,T1,'zoh');

[b1_mod,b0_mod,p_mod] = zoh_coeffs(T1,a,K);
fprintf('=== 4.1 Recovered plant ===\n');
fprintf('a = %.4f rad/s,   K = %.4f\n', a, K);
disp(Gz_check)
fprintf('Given   G(z): b1 = %.4f, b0 = %.4f (zero at z = -%.4f), pole %.4f\n', b1_given,b0_given,zero_given,pz_given);
fprintf('Model   G(z): b1 = %.4f, b0 = %.4f (zero at z = -%.4f), pole %.4f\n', b1_mod,b0_mod,b0_mod/b1_mod,p_mod);
fprintf('Zero mismatch: given %.4f vs model %.4f  (%.1f%%).  Using the %s G(z) for the T=1 s analyses.\n', ...
    zero_given, b0_mod/b1_mod, 100*(zero_given-b0_mod/b1_mod)/(b0_mod/b1_mod), ternary(USE_GIVEN_AT_T1,'GIVEN','MODEL'));

if USE_GIVEN_AT_T1
    b1T = b1_given; b0T = b0_given; pT = pz_given;     % plant used for T = 1 s sections
else
    b1T = b1_mod;   b0T = b0_mod;   pT = p_mod;
end
Gz1 = tf([b1T b0T],[1 -(1+pT) pT], T1);

%% 4.2 Analog gain range  (Objective 2)
zeta_min      = -log(OS_spec/100)/sqrt(pi^2 + log(OS_spec/100)^2);
Kp_max_analog = a^2/(4*K*zeta_min^2);
Kp_min_analog = a^2/(4*K*1^2);
ts_analog     = 8/a;

fprintf('\n=== 4.2 Analog design ===\n');
fprintf('zeta_min = %.4f\n', zeta_min);
fprintf('Analog Kp range: [%.4f, %.4f],  ideal ts = %.4f s\n', Kp_min_analog, Kp_max_analog, ts_analog);

%% T = 1 s design boundaries (used by every T = 1 s section below)
[Kp1_min,Kp1_bi] = breakaway_gain(b1T,b0T,pT);
Kp1_max  = zeta_crossing_gain(b1T,b0T,pT,zeta_min,Kp1_min*1.0005,Kp1_bi*0.999);
Kp1_stab = marginal_stability_gain(b0T,pT);
Kp1_mid  = 0.5*(Kp1_max + Kp1_stab);                 % intermediate gain for Nyquist
[OS1min,ts1min] = evaluate_design(b1T,b0T,pT,Kp1_min,T1);
[OS1max,ts1max] = evaluate_design(b1T,b0T,pT,Kp1_max,T1);
fprintf('\n=== T = 1 s design boundaries (%s G(z)) ===\n', ternary(USE_GIVEN_AT_T1,'given','model'));
fprintf('Kp,min  = %.4f   (OS %.2f%%, ts %.1f s)\n', Kp1_min, OS1min, ts1min);
fprintf('Kp,max  = %.4f   (OS %.2f%%, ts %.1f s)\n', Kp1_max, OS1max, ts1max);
fprintf('Kp,stab = %.4f   (margin %.2fx over Kp,max)\n', Kp1_stab, Kp1_stab/Kp1_max);

%% 4.3-4.7 Digital design for T = 0.1:0.1:1.0  (Objectives 3-4) -> Table 1
Ts = (1:10)/10;               % exact 0.1 ... 1.0
n  = numel(Ts);
[~,iT1]  = min(abs(Ts-1.0));
[~,iT05] = min(abs(Ts-0.5));

p_col=zeros(n,1); b1_col=zeros(n,1); b0_col=zeros(n,1);
Kp_min_col=zeros(n,1); Kp_max_col=zeros(n,1); Kp_stab_col=zeros(n,1);
OSmin=zeros(n,1); tsmin=zeros(n,1); OSmax=zeros(n,1); tsmax=zeros(n,1);
for i = 1:n
    T = Ts(i);
    [b1,b0,p] = zoh_coeffs(T,a,K);
    [Kp_min,Kp_breakin] = breakaway_gain(b1,b0,p);
    Kp_max  = zeta_crossing_gain(b1,b0,p,zeta_min,Kp_min*1.0005,Kp_breakin*0.999);
    Kp_stab = marginal_stability_gain(b0,p);
    [OSmin(i),tsmin(i)] = evaluate_design(b1,b0,p,Kp_min,T);
    [OSmax(i),tsmax(i)] = evaluate_design(b1,b0,p,Kp_max,T);
    p_col(i)=p; b1_col(i)=b1; b0_col(i)=b0;
    Kp_min_col(i)=Kp_min; Kp_max_col(i)=Kp_max; Kp_stab_col(i)=Kp_stab;
end
Table1 = table(Ts(:), p_col, b1_col, b0_col, Kp_min_col, Kp_max_col, Kp_stab_col, ...
    Kp_stab_col./Kp_max_col, OSmin, tsmin, OSmax, tsmax, ...
    'VariableNames', {'T','p','b1','b0','Kp_min','Kp_max','Kp_stab','margin_ratio', ...
                      'OS_at_Kpmin','ts_at_Kpmin','OS_at_Kpmax','ts_at_Kpmax'});
fprintf('\n=== Table 1: digital design across T (recovered continuous model) ===\n');
disp(Table1)

%% 4.8 Empirical Kp(T) laws  (Objective 5)
[Kp0_max,c_max,R2_max] = fit_inverse_law(Table1.T, Table1.Kp_max);
[Kp0_min,c_min,R2_min] = fit_inverse_law(Table1.T, Table1.Kp_min);
fprintf('\n=== 4.8 Empirical fits ===\n');
fprintf('Kp,max(T) ~= %.4f/(1+%.4f T),  R^2 = %.5f\n', Kp0_max, c_max, R2_max);
fprintf('Kp,min(T) ~= %.4f/(1+%.4f T),  R^2 = %.5f\n', Kp0_min, c_min, R2_min);
theory_c = a/2;                              % ZOH half-sample-delay heuristic
r2f = @(y,yh) 1 - sum((y-yh).^2)/sum((y-mean(y)).^2);
fprintf('Delay heuristic c=a/2=%.4f:  R2(Kp,max)=%.4f   R2(Kp,min)=%.4f\n', theory_c, ...
    r2f(Table1.Kp_max, Kp_max_analog./(1+theory_c*Table1.T)), ...
    r2f(Table1.Kp_min, Kp_min_analog./(1+theory_c*Table1.T)));

%% 5.2 Figure 1: digital root locus at T = 1 s with the three design points
Kp_sweep = logspace(-3, log10(Kp1_stab*1.15), 4000);
locus = zeros(numel(Kp_sweep),2);
for i = 1:numel(Kp_sweep)
    r = roots([1, Kp_sweep(i)*b1T-(1+pT), pT+Kp_sweep(i)*b0T]);
    locus(i,:) = sort(r(:)).';
end
figure('Name','Figure 1: root locus T=1s'); hold on;
th = linspace(0,2*pi,400);
plot(cos(th),sin(th),'k--');
plot(real(locus(:,1)),imag(locus(:,1)),'b-',real(locus(:,2)),imag(locus(:,2)),'b-');
plot([1 pT],[0 0],'kx','MarkerSize',10,'LineWidth',2);
plot(-b0T/b1T,0,'ko','MarkerSize',8);
for Kp = [Kp1_min, Kp1_max, Kp1_stab]
    r = roots([1, Kp*b1T-(1+pT), pT+Kp*b0T]);
    zc = r(imag(r) >= 0); zc = zc(1);
    plot(real(zc), imag(zc),'r^','MarkerSize',9,'MarkerFaceColor','r');
    if imag(zc) > 1e-9, plot(real(zc),-imag(zc),'r^','MarkerSize',9,'MarkerFaceColor','r'); end
end
axis equal; grid on; xlim([-1.7 1.3]); ylim([-1.5 1.5]);
xlabel('Real part'); ylabel('Imaginary part');
title('Digital root locus, T = 1 s: breakaway, 5% OS and marginal-stability points');

%% 5.3 Figures 2-3: Kp boundaries and settling time vs T
figure('Name','Figure 2: Kp vs T');
Tfine = linspace(0.08,1.05,200);
semilogy(Table1.T,Table1.Kp_stab,'d-','Color',[0.58 0.4 0.74]); hold on;
semilogy(Table1.T,Table1.Kp_max,'r^-');
semilogy(Table1.T,Table1.Kp_min,'s-','Color',[0.17 0.63 0.17]);
semilogy(Tfine,Kp0_max./(1+c_max*Tfine),'r--');
semilogy(Tfine,Kp0_min./(1+c_min*Tfine),'--','Color',[0.17 0.63 0.17]);
yline(Kp_max_analog,':r'); yline(Kp_min_analog,':g');
grid on; xlabel('Sampling time T (s)'); ylabel('Proportional gain Kp (log scale)');
legend('Kp,stability','Kp,max (\zeta=0.69)','Kp,min (\zeta=1)','fit max','fit min','Location','best');
title('Design boundaries and the true stability limit vs sampling time');

figure('Name','Figure 3: ts vs T');
plot(Table1.T,Table1.ts_at_Kpmax,'r^-'); hold on;
plot(Table1.T,Table1.ts_at_Kpmin,'s-','Color',[0.17 0.63 0.17]);
yline(ts_spec,'k-','LineWidth',1.2);
grid on; xlabel('Sampling time T (s)'); ylabel('2% settling time t_s (s)');
legend('t_s @ Kp,max (\zeta=0.69)','t_s @ Kp,min (\zeta=1)','Spec','Location','northwest');
title('Simulated settling time vs sampling time');

%% 6.1 Distance to instability
[mn,imn] = min(Table1.margin_ratio); [mx,imx] = max(Table1.margin_ratio);
fprintf('\n=== 6.1 Gain margin to instability (Kp,stab/Kp,max) ===\n');
fprintf('min %.2fx at T=%.1f,  max %.2fx at T=%.1f\n', mn, Table1.T(imn), mx, Table1.T(imx));

%% 6.2 Critical-damping settling-time correction (Table 2)
ts_naive = zeros(n,1); ts_corr = zeros(n,1);
for i = 1:n
    [b1,b0,p] = zoh_coeffs(Table1.T(i),a,K);
    r = roots([1, Table1.Kp_min(i)*b1-(1+p), p+Table1.Kp_min(i)*b0]);
    sigma = -log(abs(r(1)))/Table1.T(i);
    ts_naive(i) = 4/sigma;  ts_corr(i) = 5.8339/sigma;   % repeated real pole: 5.83/sigma
end
Table2 = table(Table1.T, ts_naive, ts_corr, Table1.ts_at_Kpmin, ...
    'VariableNames',{'T','ts_naive_4_over_sigma','ts_corrected_5p83_over_sigma','ts_simulated'});
fprintf('\n=== Table 2: settling time at the Kp,min boundary ===\n'); disp(Table2)

%% 6.4 Figure 4: step-response waterfall at T = 1 s (data kept for 6.7-6.8)
Kp_vals_w = logspace(log10(0.1*Kp1_min), log10(Kp1_stab*1.35), 70);
figure('Name','Figure 4: waterfall'); hold on;
wf_y = cell(numel(Kp_vals_w),1); wf_t = cell(numel(Kp_vals_w),1);
for k = 1:numel(Kp_vals_w)
    [y,t] = step(feedback(Kp_vals_w(k)*Gz1,1), 22);
    wf_y{k} = y; wf_t{k} = t;
    plot3(t, log10(Kp_vals_w(k))*ones(size(t)), max(min(y,3),-3), 'Color',[0.2 0.2 0.6]);
end
grid on; view(-60,22);
xlabel('Time (s)'); ylabel('log_{10}(K_p)'); zlabel('y(t)');
title('Digital step-response waterfall, T = 1 s');

%% 6.5 Figure 5: overshoot heatmap %OS(Kp,T)  (recovered continuous model)
Kp_grid = logspace(log10(0.15), log10(12), 200);
T_grid  = linspace(0.05, 1.2, 180);
OS_grid = zeros(numel(T_grid), numel(Kp_grid));
for i = 1:numel(T_grid)
    [b1,b0,p] = zoh_coeffs(T_grid(i),a,K);
    for j = 1:numel(Kp_grid)
        r = roots([1, Kp_grid(j)*b1-(1+p), p+Kp_grid(j)*b0]);
        if all(abs(imag(r)) < 1e-9)
            OS_grid(i,j) = 0;                       % real poles: no overshoot
        else
            zc = r(imag(r) > 0); zc = zc(1);
            if abs(zc) >= 1
                OS_grid(i,j) = 100;                 % unstable
            else
                ze = effective_zeta(zc);
                if ze <= 0, OS_grid(i,j) = 100;
                else, OS_grid(i,j) = min(100, 100*exp(-ze*pi/sqrt(max(1e-9,1-ze^2)))); end
            end
        end
    end
end
figure('Name','Figure 5: %OS heatmap');
pcolor(Kp_grid,T_grid,OS_grid); shading flat; hold on; set(gca,'XScale','log');
colormap(flipud(hot)); cb = colorbar; cb.Label.String = 'Percent overshoot (%)';
contour(Kp_grid,T_grid,OS_grid,[0.5 0.5],'LineColor','g','LineWidth',1.5);
contour(Kp_grid,T_grid,OS_grid,[5 5],'LineColor','w','LineWidth',1.5);
plot(Table1.Kp_min,Table1.T,'g.','MarkerSize',10);
plot(Table1.Kp_max,Table1.T,'w.','MarkerSize',10);
xlabel('K_p (log scale)'); ylabel('Sampling time T (s)');
title('Overshoot sensitivity, %OS(K_p,T)');

%% 6.7 Mine the waterfall: continuous %OS(Kp) and ts(Kp) at T = 1 s
nW = numel(Kp_vals_w);
OS_w = nan(nW,1); ts_w = nan(nW,1); unstable_w = Kp_vals_w(:) > Kp1_stab;
for k = 1:nW
    if unstable_w(k), continue, end
    y = wf_y{k}; t = wf_t{k};
    OS_w(k) = max(0,(max(y)-1)*100);
    outside = find(abs(y-1) > 0.02);
    if isempty(outside), ts_w(k) = 0; else, ts_w(k) = t(min(outside(end)+1,numel(t))); end
end
[ts_min_w,imin_w] = min(ts_w);
fprintf('\n=== 6.7 Waterfall-mined ts(Kp) ===\n');
fprintf('Minimum ts = %.1f s at Kp = %.3f (OS there = %.2f%%)\n', ts_min_w, Kp_vals_w(imin_w), OS_w(imin_w));
figure('Name','Figure 6: OS(Kp), ts(Kp)');
yyaxis left;  semilogx(Kp_vals_w(~unstable_w),OS_w(~unstable_w),'r.-'); ylabel('Percent overshoot (%)');
yyaxis right; semilogx(Kp_vals_w(~unstable_w),ts_w(~unstable_w),'b.-'); ylabel('2% settling time (s)'); hold on;
plot(Kp_vals_w(imin_w),ts_min_w,'p','MarkerSize',14,'MarkerFaceColor',[0.95 0.8 0.1]);
xlabel('K_p (log scale)'); grid on;
title('%OS(K_p) and t_s(K_p) at T = 1 s, mined from the waterfall');

%% 6.8 Validate the waterfall against the pole prediction
fprintf('\n=== 6.8 Empirical vs theoretical decay rate / frequency ===\n');
for Kc = [0.7*Kp1_min, Kp1_max, Kp1_mid, 0.93*Kp1_stab]
    [~,k] = min(abs(Kp_vals_w-Kc));
    Kp = Kp_vals_w(k); y = wf_y{k}; t = wf_t{k}; dev = y-1;
    r = roots([1, Kp*b1T-(1+pT), pT+Kp*b0T]);
    if all(abs(imag(r)) < 1e-9)
        fprintf('Kp=%.3f: poles real (below breakaway) -> no oscillation to measure.\n',Kp); continue
    end
    zc = r(imag(r)>0); zc = zc(1);  s = log(zc)/T1;
    ex = find_local_extrema(dev);
    sig_emp = NaN; wd_emp = NaN;
    if numel(ex) >= 2
        am = abs(dev(ex)); ra = am(1:end-1)./am(2:end); ra = ra(isfinite(ra) & ra>1);
        if ~isempty(ra), sig_emp = mean(log(ra))/T1; end
        pk = ex(dev(ex)>0);
        if numel(pk) >= 2, wd_emp = 2*pi/mean(diff(t(pk))); end
    end
    fprintf('Kp=%.3f: sigma theory %.4f / emp %.4f   wd theory %.4f / emp %.4f\n', ...
        Kp, -real(s), sig_emp, imag(s), wd_emp);
end

%% 6.9 Companion heatmap: ts(Kp,T)  (full simulation at every grid point)
Kp_g9 = logspace(log10(0.2), log10(10), 70);
T_g9  = linspace(0.08, 1.2, 60);
ts_g9 = nan(numel(T_g9), numel(Kp_g9));
for i = 1:numel(T_g9)
    [b1,b0,p] = zoh_coeffs(T_g9(i),a,K);
    Gz = tf([b1 b0],[1 -(1+p) p], T_g9(i));
    for j = 1:numel(Kp_g9)
        try
            inf9 = stepinfo(feedback(Kp_g9(j)*Gz,1),'SettlingTimeThreshold',0.02);
            ts_g9(i,j) = min(inf9.SettlingTime, 20);
        catch
            ts_g9(i,j) = NaN;
        end
    end
end
plot9 = ts_g9; plot9(isnan(plot9)) = 20;      % never settles -> cap at 20 s
figure('Name','Figure 7: ts(Kp,T) heatmap');
pcolor(Kp_g9,T_g9,plot9); shading flat; hold on; set(gca,'XScale','log');
colormap(flipud(bone)); cb = colorbar; cb.Label.String = '2% settling time (s), capped at 20 s';
contour(Kp_g9,T_g9,plot9,[5 7 9 12],'LineColor','w','ShowText','on');
xlabel('K_p (log scale)'); ylabel('Sampling time T (s)');
title('Settling-time landscape t_s(K_p,T)');

%% 6.10 Overshoot sensitivity field d(%OS)/dKp (reuses the Sec. 6.5 grid)
dOS = zeros(size(OS_grid));
for i = 1:size(OS_grid,1), dOS(i,:) = gradient(OS_grid(i,:), Kp_grid); end
figure('Name','Figure 8: sensitivity field');
pcolor(Kp_grid,T_grid,min(max(dOS,0),30)); shading flat; hold on;
set(gca,'XScale','log'); xlim([0.3 6]); colormap(parula);
cb = colorbar; cb.Label.String = 'd(%OS)/dK_p  (clipped at 30)';
contour(Kp_grid,T_grid,min(max(dOS,0),30),[5 10 15 20 25],'LineColor','w','ShowText','on');
ridge = zeros(size(T_grid));
for i = 1:numel(T_grid), [~,j] = max(dOS(i,:)); ridge(i) = Kp_grid(j); end
plot(ridge,T_grid,'r--','LineWidth',1.6);
xlabel('K_p (log scale)'); ylabel('Sampling time T (s)');
title('Where is overshoot most gain-sensitive?  d(%OS)/dK_p');

%% ====================== PART II: EXTENDED ANALYSIS ======================

%% 7.1 Graphical zgrid workflow, T = 1 s
% Interactive part: click the locus where it meets the dashed zeta=0.69
% curve; the datatip shows Gain, Damping, Overshoot.
h = rlocusplot(Gz1);
hold on; zgrid(0.69,[]);
opts = getoptions(h);
opts.Title.String  = 'Root locus with \zeta = 0.69 contour, T = 1 s';
opts.XLabel.String = 'Real axis';  opts.YLabel.String = 'Imaginary axis';
setoptions(h,opts);

[rl,kl] = rlocus(Gz1);                     % programmatic cross-check of the datatip
zeta_vec = nan(1,numel(kl));
for q = 1:numel(kl)
    zc = rl(imag(rl(:,q)) > 1e-9, q);
    if ~isempty(zc), zeta_vec(q) = effective_zeta(zc(1)); end
end
[~,iq] = min(abs(zeta_vec - zeta_min));
fprintf('\n=== 7.1 Graphical vs analytic gain, T = 1 s ===\n');
fprintf('Kp,max  graphical %.4f | analytic %.4f | diff %.2f%%\n', kl(iq), Kp1_max, 100*(kl(iq)-Kp1_max)/Kp1_max);
fprintf('Kp,min  analytic  %.4f  (click the zeta=1 breakaway point by hand to read it graphically)\n', Kp1_min);

%% 7.2 Control System Designer walkthrough, T = 0.5 s
[b1c,b0c,pc] = zoh_coeffs(0.5,a,K);
Gz_csd = tf([b1c b0c],[1 -(1+pc) pc], 0.5);
% controlSystemDesigner('rlocus', Gz_csd)    % <- uncomment to open the app
% Steps: (1) app opens with the plant's poles/zero
%   (2) right-click plot > Grid > enable damping grid
%   (3) right-click > Design Requirements > New > Damping Ratio = 0.69
%   (4) drag the closed-loop pole square onto the zeta = 0.69 line
%   (5) read the gain in the Compensator Editor
%   (6) drag to the real-axis breakaway point and read that gain
fprintf('\n=== 7.2 Targets for T = 0.5 s ===\n');
fprintf('Plant poles z = 1, %.4f ; zero z = %.4f\n', pc, -b0c/b1c);
fprintf('Expect Kp ~ %.4f at the zeta=1 breakaway, Kp ~ %.4f at zeta=0.69\n', Table1.Kp_min(iT05), Table1.Kp_max(iT05));

%% 7.3 Gain vs sampling time (Q3), from the fitted laws
Kp_max_fit = Kp0_max./(1+c_max*Table1.T);
Kp_min_fit = Kp0_min./(1+c_min*Table1.T);
Table9_3 = table(Table1.T, Kp_max_fit-Kp_min_fit, Kp_max_fit/Kp_max_analog*100, ...
    Kp_min_fit/Kp_min_analog*100, Table1.Kp_max-Table1.Kp_min, ...
    'VariableNames',{'T','width_from_fits','pct_Kpmax_of_analog','pct_Kpmin_of_analog','width_from_Table1'});
fprintf('\n=== 7.3 Admissible interval vs T ===\n'); disp(Table9_3)
figure('Name','Figure 9: interval width and % of analog');
subplot(1,2,1); plot(Table1.T,Table9_3.width_from_fits,'o-'); grid on; xlabel('T (s)'); ylabel('K_{p,max} - K_{p,min}');
subplot(1,2,2); plot(Table1.T,Table9_3.pct_Kpmax_of_analog,'r^-',Table1.T,Table9_3.pct_Kpmin_of_analog,'s-g'); grid on;
yline(100,':'); xlabel('T (s)'); ylabel('% of analog gain'); legend('K_{p,max}','K_{p,min}');

%% 8.1 Gain margin, phase margin, crossover frequency vs T
% [Gm,Pm,Wcg,Wcp] = margin(L): Wcg = frequency of the gain margin (phase = -180 deg),
% Wcp = frequency of the phase margin (|L| = 1), i.e. the GAIN crossover wgc.
PM_max=zeros(n,1); GMdB_max=zeros(n,1); GM_max=zeros(n,1); wgc_max=zeros(n,1);
PM_min=zeros(n,1); GMdB_min=zeros(n,1); GM_min=zeros(n,1); wgc_min=zeros(n,1);
for i = 1:n
    [b1,b0,p] = zoh_coeffs(Table1.T(i),a,K);
    Gz = tf([b1 b0],[1 -(1+p) p], Table1.T(i));
    [GM,PM,~,Wcp] = margin(Table1.Kp_max(i)*Gz);
    GM_max(i)=GM; GMdB_max(i)=20*log10(GM); PM_max(i)=PM; wgc_max(i)=Wcp;
    [GM,PM,~,Wcp] = margin(Table1.Kp_min(i)*Gz);
    GM_min(i)=GM; GMdB_min(i)=20*log10(GM); PM_min(i)=PM; wgc_min(i)=Wcp;
end
fprintf('\n=== 8.1 Margins at Kp,max ===\n');
disp(table(Table1.T,Table1.Kp_max,PM_max,GMdB_max,GM_max,wgc_max,'VariableNames',{'T','Kp','PM_deg','GM_dB','GM_lin','wgc'}))
fprintf('=== 8.1 Margins at Kp,min ===\n');
disp(table(Table1.T,Table1.Kp_min,PM_min,GMdB_min,GM_min,wgc_min,'VariableNames',{'T','Kp','PM_deg','GM_dB','GM_lin','wgc'}))
fprintf('Check: GM_lin x Kp,max at T=1 = %.4f   (Kp,stab = %.4f)\n', GM_max(iT1)*Table1.Kp_max(iT1), Table1.Kp_stab(iT1));
figure('Name','Figure 10: margins vs T');
yyaxis left;  plot(Table1.T,PM_min,'s-g',Table1.T,PM_max,'^-r'); ylabel('Phase margin (deg)'); ylim([0 90]);
yyaxis right; semilogy(Table1.T,GM_max,'d--'); ylabel('Gain margin, linear (Kp,max)');
xlabel('Sampling time T (s)'); grid on; legend('PM @ Kp,min','PM @ Kp,max','GM @ Kp,max','Location','best');
title('Phase margin and gain margin vs T');

%% 8.2 Nyquist at T = 1 s: design, intermediate and marginal gains
Kps_nyq = [Kp1_max, Kp1_mid, Kp1_stab];
figure('Name','Figure 11: Nyquist');
nyquistplot(Kps_nyq(1)*Gz1, Kps_nyq(2)*Gz1, Kps_nyq(3)*Gz1);   % auto-indents around z = 1
legend(sprintf('K_p = %.4f (design)',Kps_nyq(1)), sprintf('K_p = %.4f (intermediate)',Kps_nyq(2)), ...
       sprintf('K_p = %.4f (marginal)',Kps_nyq(3)),'Location','best');
title('Nyquist diagrams of L(z), T = 1 s');

% Open loop has no poles outside the unit circle (P = 0), so clockwise
% encirclements of -1 equal the number of unstable closed-loop poles.
fprintf('\n=== 8.2 Encirclements of -1, T = 1 s ===\n');
zz = exp(1j*linspace(1e-4,pi,6000));
for Kp = Kps_nyq
    L  = Kp*(b1T*zz+b0T)./((zz-1).*(zz-pT));
    rc = roots([1, Kp*b1T-(1+pT), pT+Kp*b0T]);
    N  = sum(abs(rc) > 1+1e-6);
    on = any(abs(abs(rc)-1) < 1e-4);
    fprintf('Kp=%.4f: closest approach to -1 = %.5f | N = %d clockwise encirclements%s\n', ...
        Kp, min(abs(L+1)), N, ternary(on,'  (locus passes through -1: marginal)',''));
end

%% 9.1 Closed-loop bandwidth vs Nyquist frequency
bw = zeros(n,1); wnyq = pi./Table1.T(:);
for i = 1:n
    [b1,b0,p] = zoh_coeffs(Table1.T(i),a,K);
    Gz = tf([b1 b0],[1 -(1+p) p], Table1.T(i));
    bw(i) = bandwidth(feedback(Table1.Kp_max(i)*Gz,1));
end
Table11 = table(Table1.T,Table1.Kp_max,bw,wnyq,wnyq./bw,'VariableNames',{'T','Kp_max','wb','pi_over_T','ratio'});
fprintf('\n=== 9.1 Bandwidth vs Nyquist frequency ===\n'); disp(Table11)
figure('Name','Figure 12: bandwidth vs Nyquist');
yyaxis left;  semilogy(Table1.T,bw,'o-',Table1.T,wnyq,'s-'); ylabel('Frequency (rad/s)');
yyaxis right; plot(Table1.T,wnyq./bw,'^:'); hold on; yline(10,'--'); yline(6,'--'); ylabel('Ratio (\pi/T)/\omega_b');
xlabel('Sampling time T (s)'); grid on; legend('\omega_b','\pi/T','ratio','Location','best');
title('Closed-loop bandwidth vs Nyquist frequency');

%% 9.2 ZOH zero migration  z0(T) = -b0/b1
Tg = linspace(0.05,1.2,400);
[b1g,b0g,pg] = zoh_coeffs(Tg,a,K);
figure('Name','Figure 13: zero migration');
plot(Tg,-b0g./b1g,'LineWidth',2); hold on; plot(Tg,pg,'r--'); yline(-1,':k');
plot(T1,-zero_given,'ko','MarkerFaceColor','y');            % the zero actually written on the lab sheet
xlabel('T (s)'); ylabel('z-plane location'); grid on;
legend('ZOH zero z_0(T) of recovered model','Open-loop pole p(T)','','Given zero at T=1 s');
title('ZOH zero migration vs sampling time');

%% 10.1 Pole sensitivity |dz/dKp| at T = 1 s
[s_des,z_des] = pole_sensitivity(Kp1_max, b1T,b0T,pT);
[s_mar,z_mar] = pole_sensitivity(Kp1_stab,b1T,b0T,pT);
fprintf('\n=== 10.1 Pole sensitivity ===\n');
fprintf('Kp=%.4f: z=%.4f%+.4fi  |dz/dKp|=%.5f\n', Kp1_max, real(z_des),imag(z_des),abs(s_des));
fprintf('Kp=%.4f: z=%.4f%+.4fi  |dz/dKp|=%.5f\n', Kp1_stab,real(z_mar),imag(z_mar),abs(s_mar));
fprintf('Ratio design/marginal = %.3f\n', abs(s_des)/abs(s_mar));

%% 10.2 Robustness to +-10% errors in a and K  (continuous model, nominal Kp)
Kp_nom = Table1.Kp_max(iT1);                  % nominal gain of the continuous-model design
combos = [1.1 1.1; 1.1 0.9; 0.9 1.1; 0.9 0.9];
OSr = zeros(4,1); tsr = zeros(4,1);
fprintf('\n=== 10.2 Robustness, T = 1 s, nominal Kp = %.4f ===\n', Kp_nom);
for i = 1:4
    [b1_,b0_,p_] = zoh_coeffs(1.0, a*combos(i,1), K*combos(i,2));
    [OSr(i),tsr(i)] = evaluate_design(b1_,b0_,p_,Kp_nom,1.0);
    fprintf('a%+.0f%%  K%+.0f%%:  OS = %6.2f%%   ts = %.1f s   %s\n', ...
        (combos(i,1)-1)*100,(combos(i,2)-1)*100,OSr(i),tsr(i), ...
        ternary(OSr(i)<=OS_spec+OS_tol && tsr(i)<=ts_spec,'meets spec','VIOLATES spec'));
end
fprintf('Worst case: OS = %.2f%%, ts = %.1f s\n', max(OSr), max(tsr));

%% 11. Lead compensator, T = 1 s:  C(z) = K (z - z_lead)/(z - p_lead)
% The zero cancels the plant pole (z_lead = pT); the loop then has the same
% form as the proportional design with the pole pT replaced by the free
% p_lead, so breakaway_gain / zeta_crossing_gain are reused with p_lead.
% Pick the least aggressive p_lead whose FULL simulation meets
% OS <= spec and ts <= ts_lead - lead_margin.
z_lead = pT;

% Fairness check first: how fast can a PROPORTIONAL controller already settle
% inside its own admissible window [Kp,min, Kp,max]?  (ts is quantised to
% whole sampling periods, so the best value holds over a band of Kp.)
Kscan   = linspace(Kp1_min, Kp1_max, 400);
ts_scan = nan(size(Kscan));  OS_scan = nan(size(Kscan));
for q = 1:numel(Kscan)
    [OS_scan(q),ts_scan(q)] = evaluate_design(b1T,b0T,pT,Kscan(q),T1);
end
ts_prop_best = min(ts_scan);
idx_best     = find(ts_scan == ts_prop_best);
Kp_best      = Kscan(idx_best(max(1,round(numel(idx_best)/2))));   % middle of the best band
[OS_best,~]  = evaluate_design(b1T,b0T,pT,Kp_best,T1);
fprintf('\n=== 11. Best proportional-only design ===\n');
fprintf('Fastest ts = %.1f s for Kp in [%.4f, %.4f]; mid-band Kp = %.4f gives OS = %.2f%%\n', ...
    ts_prop_best, Kscan(idx_best(1)), Kscan(idx_best(end)), Kp_best, OS_best);
fprintf('A lead compensator is only worthwhile if it settles faster than %.1f s.\n', ts_prop_best);

cand = 0.15:-0.05:-0.40;
Kc_list = nan(size(cand)); OSc = nan(size(cand)); tsc = nan(size(cand));
for q = 1:numel(cand)
    try
        [km,kb] = breakaway_gain(b1T,b0T,cand(q));
        Kc_list(q) = zeta_crossing_gain(b1T,b0T,cand(q),zeta_min,km*1.0005,kb*0.999);
        Cq = tf(Kc_list(q)*[1 -z_lead],[1 -cand(q)],T1);
        iq2 = stepinfo(feedback(Cq*Gz1,1),'SettlingTimeThreshold',0.02);
        OSc(q) = iq2.Overshoot; tsc(q) = iq2.SettlingTime;
    catch
        % bracket fails for very negative p_lead: skip that candidate
    end
end
fprintf('\n=== 11. Lead search over p_lead (z_lead = %.3f) ===\n', z_lead);
disp(table(cand(:),Kc_list(:),OSc(:),tsc(:),'VariableNames',{'p_lead','K','OS_pct','ts_s'}))
ok = find(OSc <= OS_spec+OS_tol & tsc <= ts_lead-lead_margin, 1, 'first');
if isempty(ok)
    [~,ok] = min(tsc);
    fprintf('WARNING: no candidate met ts <= %g s with %g s margin; using the fastest one.\n', ts_lead, lead_margin);
end
p_lead = cand(ok);  K_lead = Kc_list(ok);
C = tf(K_lead*[1 -z_lead],[1 -p_lead],T1);
Gcl_lead = feedback(C*Gz1,1);
info = stepinfo(Gcl_lead,'SettlingTimeThreshold',0.02);
fprintf('\nC(z) = %.4f (z - %.3f)/(z %+.3f)\n', K_lead, z_lead, -p_lead);
fprintf('Verification (full discrete step): OS = %.2f%%,  ts = %.2f s   (target ts <= %g s)\n', ...
    info.Overshoot, info.SettlingTime, ts_lead);

figure('Name','Figure 14: compensated vs proportional');
[yP,tP] = step(feedback(Kp1_max*Gz1,1),16);
[yB,tB] = step(feedback(Kp_best*Gz1,1),16);
[yC,tC] = step(Gcl_lead,16);
stairs(tP,yP,'Color',[0.5 0.5 0.5],'LineWidth',1.4); hold on;
stairs(tB,yB,'b--','LineWidth',1.4); stairs(tC,yC,'r','LineWidth',1.8);
yline(1,'k-'); yline(1.02,'k:'); yline(0.98,'k:'); xline(ts_lead,'g--'); grid on;
xlabel('Time (s)'); ylabel('Output y(t)');
legend('Proportional, K_{p,max}','Proportional, best ts','Lead compensator','Location','southeast');
title('Lead-compensated vs proportional-only step response, T = 1 s');

figure('Name','Compensated root locus');
rlocus(minreal(tf([1 -z_lead],[1 -p_lead],T1)*Gz1)); hold on; zgrid(0.69,[]);   % gain K factored out
rc = roots([1,K_lead*b1T-(1+p_lead),p_lead+K_lead*b0T]);
plot(real(rc),imag(rc),'ms','MarkerFaceColor','m');
title('Compensated root locus with \zeta = 0.69 grid (magenta = chosen closed-loop poles)');

%% 12.1 Design nomograph  (recovered continuous model, grid simulated)
Tg14  = linspace(0.1,1.0,19);
Kpg14 = linspace(0.6*min(Table1.Kp_min), 1.7*max(Table1.Kp_max), 85);
OS14 = zeros(numel(Kpg14),numel(Tg14)); ts14 = OS14;
for i = 1:numel(Tg14)
    [b1,b0,p] = zoh_coeffs(Tg14(i),a,K);
    for j = 1:numel(Kpg14)
        [OS14(j,i),ts14(j,i)] = evaluate_design(b1,b0,p,Kpg14(j),Tg14(i));
    end
end
ok14 = double((OS14 <= OS_spec) & (ts14 <= ts_spec));
figure('Name','Figure 15: nomograph');
contourf(Tg14,Kpg14,ok14,[0.5 1.5],'LineStyle','none'); hold on;
colormap(gca,[1 1 1; 0.56 0.82 0.62]); caxis([0 1]);
contour(Tg14,Kpg14,OS14,[OS_spec OS_spec],'LineColor','r','LineWidth',2);
contour(Tg14,Kpg14,ts14,[ts_spec ts_spec],'LineColor',[0.12 0.22 0.39],'LineWidth',2);
plot(Table1.T,Table1.Kp_max,'ro','MarkerFaceColor','r');
plot(Table1.T,Table1.Kp_min,'gs','MarkerFaceColor','g');
xlabel('Sampling time T (s)'); ylabel('Proportional gain K_p');
title('Design nomograph: %OS \leq 5% and t_s \leq 9 s (green = compliant)');

%% 12.2 Intersample peak, T = 1 s, nominal Kp (continuous model)
% Fine ZOH sub-stepping of the continuous plant: the control is held over
% each T, the plant is integrated exactly in between.
Ac = [0 1; 0 -a];  Bc = [0; K];  Cc = [1 0];
Nsub = 60;  Kp14 = Table1.Kp_max(iT1);
sd = c2d(ss(Ac,Bc,Cc,0), T1/Nsub, 'zoh');
x = [0;0];  t_h = 0;  y_h = 0;
for m = 0:21
    u = Kp14*(1 - Cc*x);                       % sampled error, held for one period
    for k = 1:Nsub
        x = sd.A*x + sd.B*u;
        t_h(end+1) = t_h(end) + T1/Nsub;       %#ok<SAGROW>
        y_h(end+1) = Cc*x;                     %#ok<SAGROW>
    end
end
y_s = y_h(1:Nsub:end);  t_s = t_h(1:Nsub:end);
OS_s = max(0,(max(y_s)-1)*100);  [pk,ipk] = max(y_h);  OS_c = max(0,(pk-1)*100);
fprintf('\n=== 12.2 Intersample peak, T = 1 s, Kp = %.4f ===\n', Kp14);
fprintf('OS on samples   : %.4f %%\n', OS_s);
fprintf('OS continuous   : %.4f %%  (peak at t = %.3f s)\n', OS_c, t_h(ipk));
fprintf('Gap             : %.4f percentage points\n', OS_c-OS_s);
figure('Name','Figure 16: intersample peak');
plot(t_h,y_h,'b-'); hold on; plot(t_s,y_s,'ro','MarkerFaceColor','r');
xlim([2 8]); yline(1,'k-'); grid on; xlabel('Time (s)'); ylabel('Output y(t)');
legend('Continuous (intersample) response','Sampled sequence y[n]');
title(sprintf('Intersample peak, T = 1 s, K_p = %.4f', Kp14));

fprintf('\n=== ALL SECTIONS COMPLETE ===\n');

%% ======================= Local helper functions =========================
% (Script-local functions: must stay at the end of the file.)

function idx = find_local_extrema(v)
%FIND_LOCAL_EXTREMA  Indices of local maxima/minima (no Signal Processing Toolbox).
    v = v(:).';
    s = sign(diff(v));  s(s==0) = 1;
    idx = find(diff(s) ~= 0) + 1;
end

function out = ternary(cond,a,b)
%TERNARY  Inline if/else for use inside fprintf arguments.
    if cond, out = a; else, out = b; end
end
