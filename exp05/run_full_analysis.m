function run_full_analysis()
%% run_full_analysis.m
% Master script: runs Parts 1-7, then produces every figure and CSV
% table used in the report. This file contains the plotting/export
% code and the bulk grid-sweep loops that were deliberately left out
% of the report's code listings (which show only the core numerical
% methods); everything here is a direct, mechanical use of the
% functions defined in part1-part7.
clc; close all
outFigs = fullfile(pwd,'figs'); outData = fullfile(pwd,'data');
if ~exist(outFigs,'dir'), mkdir(outFigs); end
if ~exist(outData,'dir'), mkdir(outData); end

num_GOL = [1 1.7 -1];
den_GOL = [1 0 0 -1];
closedLoopPoly = @(K) den_GOL + K*[0 num_GOL];
perturbedDen   = @(f) conv([1 -1],[1 (1+f/100) (1+f/100)^2]);
closedLoopPolyF = @(K,den_f) den_f + K*[0 num_GOL];

%% ---- Table 1: poles & zeros ----
zz = roots(num_GOL); pp = roots(den_GOL);
T1 = table([zz;pp], [repmat("Zero",numel(zz),1); repmat("Pole",numel(pp),1)], ...
    real([zz;pp]), imag([zz;pp]), abs([zz;pp]), angle([zz;pp])*180/pi, ...
    'VariableNames', {'value','type','re','im','mag','angle_deg'});
writetable(T1(:,2:end), fullfile(outData,'table1_poles_zeros.csv'))

figure('Position',[100 100 550 520]);
th = linspace(0,2*pi,400);
plot(cos(th),sin(th),'--','Color',[0.35 0.56 0.84]); hold on
scatter(real(zz),imag(zz),120,'o','MarkerEdgeColor',[0.85 0.33 0.31],'LineWidth',2)
scatter(real(pp),imag(pp),130,'x','MarkerEdgeColor',[0.91 0.54 0.05],'LineWidth',2.2)
axis equal; grid on; xlabel('Real part'); ylabel('Imaginary part')
title('Pole-Zero Map of G_{OL}(z)'); legend('Unit circle','Zeros','Poles','Location','northwest')
exportgraphics(gcf, fullfile(outFigs,'fig1_polezero.png'),'Resolution',220)

%% ---- Fig 2: root locus + stability sweep ----
Kpos = linspace(0,60,3000); Kneg = linspace(0,-60,3000);
rlocusPos = zeros(numel(Kpos),3); rlocusNeg = zeros(numel(Kneg),3);
prev = [];
for i=1:numel(Kpos)
    r = roots(closedLoopPoly(Kpos(i)));
    if ~isempty(prev), r = matchRoots(prev,r); end
    rlocusPos(i,:) = r; prev = r;
end
prev = [];
for i=1:numel(Kneg)
    r = roots(closedLoopPoly(Kneg(i)));
    if ~isempty(prev), r = matchRoots(prev,r); end
    rlocusNeg(i,:) = r; prev = r;
end
Ksweep = linspace(-15,15,601);
maxAbsRoot = arrayfun(@(K) max(abs(roots(closedLoopPoly(K)))), Ksweep);
writetable(table(Ksweep',maxAbsRoot',maxAbsRoot'<1,'VariableNames',{'K','max_abs_root','stable'}), ...
    fullfile(outData,'table2_stability_sweep.csv'))

figure('Position',[100 100 1150 520]);
subplot(1,2,1)
plot(cos(th),sin(th),'--','Color',[0.35 0.56 0.84]); hold on
plot(real(rlocusPos),imag(rlocusPos),'Color',[0.16 0.62 0.56],'LineWidth',1.6)
plot(real(rlocusNeg),imag(rlocusNeg),'Color',[0.85 0.33 0.31],'LineWidth',1.6)
scatter(real(pp),imag(pp),120,'x','MarkerEdgeColor',[0.91 0.54 0.05],'LineWidth',2)
scatter(real(zz),imag(zz),110,'o','MarkerEdgeColor','k','LineWidth',1.6)
axis equal; grid on; xlabel('Real part'); ylabel('Imaginary part'); title('Root Locus of 1+K G_{OL}(z)=0')
subplot(1,2,2)
plot(Ksweep,maxAbsRoot,'Color',[0.15 0.27 0.33],'LineWidth',1.8); hold on
yline(1,'--','Color',[0.85 0.33 0.31],'LineWidth',1.4)
ylim([0.9 8]); grid on; xlabel('Proportional gain K'); ylabel('max_i |z_i(K)|')
title('Stability Sweep')
exportgraphics(gcf, fullfile(outFigs,'fig2_rootlocus_stability.png'),'Resolution',220)

%% ---- Fig 3-6: frequency response, Ts effect, filters, sideband ----
theta = linspace(1e-3,pi-1e-4,4000); z_theta = exp(1j*theta);
Gz = polyval(num_GOL,z_theta)./polyval(den_GOL,z_theta);
magdB = 20*log10(abs(Gz)); phase = unwrap(angle(Gz))*180/pi;
writetable(table(theta',magdB',phase','VariableNames',{'theta_rad','mag_dB','phase_deg'}), ...
    fullfile(outData,'mainband_response.csv'))

figure('Position',[100 100 780 680]);
subplot(2,1,1); plot(theta,magdB,'Color',[0.11 0.21 0.34],'LineWidth',1.5); hold on
xline(2*pi/3,'-.','Color',[0.85 0.33 0.31]); xline(pi,'--','Color',[0.16 0.62 0.56])
grid on; ylabel('Magnitude (dB)'); title('Main-Band Frequency Response of G_{OL}(e^{j\theta})')
subplot(2,1,2); plot(theta,phase,'Color',[0.11 0.21 0.34],'LineWidth',1.5); hold on
xline(2*pi/3,'-.','Color',[0.85 0.33 0.31]); xline(pi,'--','Color',[0.16 0.62 0.56])
grid on; xlabel('\theta (rad)'); ylabel('Phase (deg)')
exportgraphics(gcf, fullfile(outFigs,'fig3_mainband_response.png'),'Resolution',220)

theta_c = pi/3; Ts_design = 1.0; wc = (2/Ts_design)*tan(theta_c/2);
C1_z = c2d(tf(wc,[1 wc]),Ts_design,'tustin');
C2_z = c2d(tf(wc^2,[1 sqrt(2)*wc wc^2]),Ts_design,'tustin');
[numC1,denC1] = tfdata(C1_z,'v'); [numC2,denC2] = tfdata(C2_z,'v');
H1 = polyval(numC1,z_theta)./polyval(denC1,z_theta);
H2 = polyval(numC2,z_theta)./polyval(denC2,z_theta);
casc1_dB = 20*log10(abs(H1.*Gz)); casc2_dB = 20*log10(abs(H2.*Gz));
writetable(table(theta',magdB',20*log10(abs(H1))',20*log10(abs(H2))',casc1_dB',casc2_dB', ...
    'VariableNames',{'theta_rad','plant_dB','filt1_dB','filt2_dB','casc1_dB','casc2_dB'}), ...
    fullfile(outData,'filter_comparison.csv'))

z_edge = exp(1j*(pi-1e-4));
G_edge = polyval(num_GOL,z_edge)/polyval(den_GOL,z_edge);
H1_edge = polyval(numC1,z_edge)/polyval(denC1,z_edge);
H2_edge = polyval(numC2,z_edge)/polyval(denC2,z_edge);
T3 = table(["Uncompensated G_OL(z)";"1st-order filter + G_OL(z)";"2nd-order filter + G_OL(z)"], ...
    [20*log10(abs(G_edge));20*log10(abs(H1_edge*G_edge));20*log10(abs(H2_edge*G_edge))], ...
    'VariableNames',{'configuration','mag_dB_at_boundary'});
writetable(T3, fullfile(outData,'table3_boundary_filtered.csv'))

figure('Position',[100 100 780 500]);
plot(theta,magdB,'Color',[0.15 0.27 0.33],'LineWidth',1.6); hold on
plot(theta,casc1_dB,'Color',[0.16 0.62 0.56],'LineWidth',1.6)
plot(theta,casc2_dB,'Color',[0.90 0.22 0.27],'LineWidth',1.6)
xline(pi,':','Color',[0.4 0.4 0.4])
grid on; xlabel('\theta = \omega T_s (rad)'); ylabel('Magnitude (dB)')
legend('Uncompensated |G_{OL}|','Compensated |C_1 G_{OL}|','Compensated |C_2 G_{OL}|','Band edge')
title('Cascade Response: Anti-Aliasing Filter x Plant')
exportgraphics(gcf, fullfile(outFigs,'fig5_filter_comparison.png'),'Resolution',220)

Ts_sb = 1.0; theta_full = linspace(-3*pi,3*pi,20000);
theta_full(abs(mod(theta_full+pi,2*pi)-pi)<1e-3) = [];
z_full = exp(1j*theta_full); omega_full = theta_full/Ts_sb;
Gz_full = polyval(num_GOL,z_full)./polyval(den_GOL,z_full);
H2_full = polyval(numC2,z_full)./polyval(denC2,z_full);
writetable(table(omega_full',20*log10(abs(Gz_full))',20*log10(abs(H2_full.*Gz_full))', ...
    'VariableNames',{'omega','uncompensated_dB','compensated_2nd_dB'}), ...
    fullfile(outData,'sideband_response.csv'))

%% ---- K_max_stable(f) numerical scan (Part 4/5 support) ----
fScan = linspace(-30,30,25); KmaxStable = nan(size(fScan));
for i = 1:numel(fScan)
    denf = perturbedDen(fScan(i));
    Ks = logspace(-3,log10(500),4000);
    stab = arrayfun(@(K) max(abs(roots(closedLoopPolyF(K,denf))))<1, Ks);
    idx = find(~stab,1,'first');
    if ~isempty(idx) && idx>1, KmaxStable(i) = Ks(idx-1); end
end
writetable(table(fScan',KmaxStable','VariableNames',{'f_percent','K_max_stable'}), ...
    fullfile(outData,'table_kmax_vs_f.csv'))

%% ---- PO(K,f) exact vs approximate grid (Part 4) ----
Kgrid = logspace(log10(1e-4),log10(1.0),90); fgrid = linspace(-30,30,61);
PO_exact = zeros(numel(fgrid),numel(Kgrid)); PO_approx = PO_exact;
for fi = 1:numel(fgrid)
    for ki = 1:numel(Kgrid)
        f = fgrid(fi); K = Kgrid(ki);
        denf = perturbedDen(f);
        r = roots(closedLoopPolyF(K,denf));
        [~,idx] = max(abs(r)); zd = r(idx);
        y = stepResponseClosedLoop(K,f,150);
        PO_exact(fi,ki) = percentOvershoot(y,1.0,100);
        PO_approx(fi,ki) = min(poApproxFromPole(zd),100);
    end
end
save(fullfile(outData,'po_heatmap.mat'),'Kgrid','fgrid','PO_exact','PO_approx')

figure('Position',[100 100 1500 480]);
subplot(1,3,1); pcolor(Kgrid,fgrid,PO_exact); shading interp; set(gca,'XScale','log')
colormap(gca,'parula'); colorbar; xlabel('K (log)'); ylabel('f (%)'); title('Exact PO(K,f)')
subplot(1,3,2); pcolor(Kgrid,fgrid,PO_approx); shading interp; set(gca,'XScale','log')
colorbar; xlabel('K (log)'); title('Approximate PO(K,f)')
subplot(1,3,3); pcolor(Kgrid,fgrid,abs(PO_exact-PO_approx)); shading interp; set(gca,'XScale','log')
colorbar; xlabel('K (log)'); title('|Error|')
exportgraphics(gcf, fullfile(outFigs,'fig8_po_heatmap.png'),'Resolution',220)

%% ---- Waterfall step response (Part 4) ----
[~,ix] = max(KmaxStable); f_wf = fScan(ix);
Kwf = logspace(log10(0.01),log10(1.0),40); n = 55;
Ywf = zeros(numel(Kwf),n);
for i = 1:numel(Kwf)
    Ywf(i,:) = stepResponseClosedLoop(Kwf(i),f_wf,n);
    Ywf(i,:) = max(min(Ywf(i,:),6),-6);
end
save(fullfile(outData,'waterfall.mat'),'Kwf','Ywf','f_wf')

figure('Position',[100 100 800 650]);
hold on
cmap = parula(numel(Kwf));
for i = 1:numel(Kwf)
    plot3(1:n, repmat(log10(Kwf(i)),1,n), Ywf(i,:), 'Color', cmap(i,:), 'LineWidth', 1)
end
grid on; view(-60,22); xlabel('Time step k'); ylabel('log_{10}(K)'); zlabel('y[k]')
title(sprintf('Step-Response Waterfall, f=%.0f%%, K=%.3f\\to%.2f', f_wf, min(Kwf), max(Kwf)))
exportgraphics(gcf, fullfile(outFigs,'fig10_waterfall.png'),'Resolution',220)

%% ---- Part 6: exact (alpha,K) stability region + optimal triple pole ----
alphaGrid = linspace(-2.999,-0.001,400);
Klower = -(3+alphaGrid)/2;
writetable(table(alphaGrid',Klower',zeros(size(alphaGrid))', ...
    'VariableNames',{'alpha','K_lower_exact','K_upper_exact'}), ...
    fullfile(outData,'table6_alpha_region_boundary.csv'))

phi = (1+sqrt(5))/2; r_star = phi^(1/3)-phi^(-1/3);
alpha_star = -r_star; K_star = r_star^3-1;
rootsOpt = roots([1 K_star alpha_star*K_star -(K_star+1)]);
writetable(table(r_star,alpha_star,K_star,max(abs(rootsOpt)), ...
    'VariableNames',{'r_star','alpha_star','K_star','max_abs_root'}), ...
    fullfile(outData,'table6_alpha_optimal_point.csv'))

figure('Position',[100 100 550 500]);
fill([alphaGrid fliplr(alphaGrid)],[Klower zeros(size(Klower))],[0.16 0.62 0.56], ...
    'FaceAlpha',0.28,'EdgeColor','none'); hold on
plot(alphaGrid,Klower,'Color',[0.11 0.21 0.34],'LineWidth',2)
plot(alpha_star,K_star,'p','MarkerSize',16,'MarkerFaceColor',[0.90 0.22 0.27],'MarkerEdgeColor','k')
xline(1.7,'--','Color',[0.90 0.45 0.32],'LineWidth',1.8)
grid on; xlabel('\alpha'); ylabel('K'); title('Exact (\alpha,K) Stability Region')
exportgraphics(gcf, fullfile(outFigs,'fig11_alpha_region.png'),'Resolution',220)

%% ---- Part 7: zeta(K,f) / ringing map + general pole sensitivity ----
ZETA = nan(numel(fgrid),numel(Kgrid)); RINGING = zeros(numel(fgrid),numel(Kgrid));
for fi = 1:numel(fgrid)
    for ki = 1:numel(Kgrid)
        f = fgrid(fi); K = Kgrid(ki);
        r = roots(closedLoopPolyF(K,perturbedDen(f)));
        [~,idx] = max(abs(r)); zd = r(idx);
        stable = max(abs(r))<1;
        if ~stable
            RINGING(fi,ki) = -1;
        elseif abs(imag(zd))>1e-9
            RINGING(fi,ki) = 1;
            lr = log(abs(zd)); ZETA(fi,ki) = -lr/sqrt(lr^2+angle(zd)^2);
        else
            RINGING(fi,ki) = 0; ZETA(fi,ki) = 1;
        end
    end
end
save(fullfile(outData,'zeta_ringing.mat'),'Kgrid','fgrid','ZETA','RINGING')

figure('Position',[100 100 750 550]);
pcolor(Kgrid,fgrid,ZETA); shading interp; set(gca,'XScale','log'); colorbar
xlabel('K (log)'); ylabel('f (%)'); title('\zeta(K,f)')
exportgraphics(gcf, fullfile(outFigs,'fig12_zeta_heatmap.png'),'Resolution',220)

%% ---- Part 8: anti-aliasing filter generalisation (Section 10) ----
thetaPlot = linspace(1e-3,pi-1e-6,600);
thetaC = pi/3; theta60fn = @(tc,n) 2*atan(tan(tc/2)*(1e6-1)^(1/(2*n)));
magDbFn = @(th,tc,n) -10*log10(1+(tan(th/2)./tan(tc/2)).^(2*n));

tcGrid = linspace(deg2rad(20),deg2rad(110),91);
rows = []; thetaRes = 2*pi/3;
for n = 1:4
    for tc = tcGrid
        guardDeg = rad2deg(pi-theta60fn(tc,n));
        resCostDB = magDbFn(thetaRes,tc,n);
        rows = [rows; n, rad2deg(tc), guardDeg, resCostDB]; %#ok<AGROW>
    end
end
writetable(array2table(rows,'VariableNames',{'n','theta_c_deg','guard_band_deg','resonance_cost_dB'}), ...
    fullfile(outData,'table12_filter_tradeoff_sweep.csv'))

omega_c_phys = 1.0; n_fixed = 2; rowsTs = [];
for Ts = [1.0 0.5 0.2 0.1]
    tcReq = omega_c_phys*Ts;
    guardDeg = rad2deg(pi-theta60fn(tcReq,n_fixed));
    resCostDB = magDbFn(thetaRes,tcReq,n_fixed);
    rowsTs = [rowsTs; Ts, rad2deg(tcReq), guardDeg, resCostDB]; %#ok<AGROW>
end
writetable(array2table(rowsTs,'VariableNames',{'Ts','theta_c_deg','guard_band_deg','resonance_cost_dB'}), ...
    fullfile(outData,'table13_filter_ts_sensitivity.csv'))

figure('Position',[100 100 780 520]);
cmap8 = [0.91 0.44 0.32; 0.16 0.62 0.56; 0.15 0.27 0.33; 0.45 0.04 0.70];
hold on
for n = 1:4
    plot(rad2deg(thetaPlot), max(magDbFn(thetaPlot,thetaC,n),-160), 'Color', cmap8(n,:), 'LineWidth', 1.8)
end
xline(180,'k-','LineWidth',1); grid on; ylim([-160 8])
xlabel('\theta (deg)'); ylabel('|C_n(e^{j\theta})| (dB)')
legend('n=1','n=2','n=3','n=4','Location','southwest')
title('Exact Roll-off at Fixed Corner \theta_c=60 deg (Section 10)')
exportgraphics(gcf, fullfile(outFigs,'fig16_filter_tradeoff.png'),'Resolution',220)
disp('Part 8 (filter generalisation) numeric core cross-checked against part8_filter_generalisation.m in isolation.')

disp('All figures and CSV/MAT data files written to figs/ and data/.')
end

%% ---- local functions ----
function r2 = matchRoots(prev,r)
    used = false(1,numel(r)); order = zeros(1,numel(r));
    for a = 1:numel(prev)
        d = abs(prev(a)-r); d(used) = inf;
        [~,j] = min(d); used(j) = true; order(a) = j;
    end
    r2 = r(order);
end

function y = stepResponseClosedLoop(K, f_percent, n)
    num_GOL = [1 1.7 -1];
    rho = 1 + f_percent/100;
    den_f = conv([1 -1],[1 rho rho^2]);
    a = den_f + K*[0 num_GOL]; b = K*num_GOL; a = a/a(1); bpad = [0 b];
    clipVal = 1e6; y = zeros(1,n); u = ones(1,n);
    for k = 1:n
        acc = 0;
        for j = 2:4
            if k-j+1 >= 1, acc = acc + bpad(j)*u(k-j+1) - a(j)*y(k-j+1); end
        end
        yk = acc/a(1); if ~isfinite(yk), yk = clipVal; end
        y(k) = max(min(yk,clipVal),-clipVal);
    end
end

function po = percentOvershoot(y, finalValue, cap)
    y(~isfinite(y)) = cap*10;
    po = 100*(max(y)-finalValue)/finalValue;
    po = min(max(po,0),cap);
end

function po = poApproxFromPole(zd)
    r = abs(zd); phi = abs(angle(zd));
    if phi < 1e-9 || r <= 0, po = 0; return; end
    lr = log(r);
    po = 100*exp(pi*lr/phi);
    po = min(max(po,0),200);
end
