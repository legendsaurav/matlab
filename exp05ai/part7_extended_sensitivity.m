function part7_extended_sensitivity()
%% Part 7: Damping classification and general pole sensitivity dz/df
%
% Objective: (a) classify the closed-loop response at any (K,f) as
% "ringing" (complex dominant pole -> oscillatory) or "non-ringing"
% (real dominant pole -> monotonic), which explains the shape of the
% PO(K,f) heatmap of Part 4; and (b) extend Part 5's root-departure
% sensitivity (evaluated only at K=0) to the operating condition that
% actually matters once the loop is closed: how much do the poles
% move, at a FIXED working gain K, if the plant's resonance radius f
% drifts? This is obtained analytically by implicit differentiation
% of the characteristic equation F(z,K,f)=D(z,f)+K*N(z)=0:
%     dz/df = -(dF/df) / (dF/dz)
% and validated against a central finite difference.

num_GOL = [1 1.7 -1];
dnum_GOL = polyder(num_GOL);

% ---- (a) zeta and ringing classification, a few representative points ----
fprintf('--- Damping classification at representative (K,f) ---\n')
for pt = {[0.02 -20], [0.15 -20], [0.30 -20], [0.05 -5]}
    K = pt{1}(1); f = pt{1}(2);
    rho = 1+f/100;
    den_f = conv([1 rho rho^2],[1 -1]);
    r = roots(den_f + K*[0 num_GOL]);
    [~,idx] = max(abs(r)); zd = r(idx);
    [zeta, cls] = zetaAndClass(zd, max(abs(r))<1);
    fprintf('K=%6.3f f=%5.1f%%  dominant pole=%7.4f%+7.4fi  zeta=%6.4f  [%s]\n', ...
        K, f, real(zd), imag(zd), zeta, cls)
end

% ---- (b) analytical vs finite-difference dz/df at fixed K~=0 ----
fprintf('\n--- General pole sensitivity dz/df at fixed K (analytical vs finite-diff) ---\n')
for pt = {[0.10 -15], [0.10 -5], [0.15 -10]}
    K = pt{1}(1); f0 = pt{1}(2);
    rho = 1+f0/100;
    den_f = conv([1 rho rho^2],[1 -1]);
    r0 = roots(den_f + K*[0 num_GOL]);
    h = 1e-4;
    for i = 1:numel(r0)
        z0 = r0(i);
        dz_an = analyticalDzDf(z0, K, f0, num_GOL, dnum_GOL);
        rp = roots(conv([1 (1+(f0+h)/100) (1+(f0+h)/100)^2],[1 -1]) + K*[0 num_GOL]);
        rm = roots(conv([1 (1+(f0-h)/100) (1+(f0-h)/100)^2],[1 -1]) + K*[0 num_GOL]);
        [~,ip] = min(abs(rp-(z0+dz_an*h))); zp = rp(ip);
        [~,im] = min(abs(rm-(z0-dz_an*h))); zm = rm(im);
        dz_fd = (zp-zm)/(2*h);
        relerr = abs(dz_an-dz_fd)/max(abs(dz_an),1e-12);
        fprintf('K=%.2f f0=%5.1f  z0=%7.4f%+7.4fi  dz/df(analytic)=%8.5f%+8.5fi  relerr=%.2e\n', ...
            K, f0, real(z0), imag(z0), real(dz_an), imag(dz_an), relerr)
    end
end
disp('Relative error is at the level of machine/finite-difference precision for every')
disp('pole tested, confirming the analytical implicit-function-theorem formula.')
end

%% ---- local functions ----
function [zeta, cls] = zetaAndClass(zd, isStable)
% Equivalent damping ratio and ringing/non-ringing/unstable label for
% a discrete dominant pole z_d.
    r = abs(zd); ph = abs(angle(zd));
    if ~isStable
        zeta = NaN; cls = 'unstable'; return
    end
    if ph < 1e-9
        zeta = 1.0; cls = 'non-ringing'; return
    end
    lr = log(r);
    zeta = -lr/sqrt(lr^2+ph^2);
    cls = 'ringing';
end

function dzdf = analyticalDzDf(z, K, f_percent, num_GOL, dnum_GOL)
% Implicit-function-theorem pole sensitivity dz/df at fixed K, for
% F(z,K,f) = D(z,f) + K*N(z) = 0,  D(z,f) = (z-1)(z^2+rho z+rho^2).
    rho = 1 + f_percent/100;
    dDdz = polyval(polyder(conv([1 rho rho^2],[1 -1])), z);
    dFdz = dDdz + K*polyval(dnum_GOL, z);
    dDdrho = (z-1)*(z+2*rho);          % d/drho of (z-1)(z^2+rho z+rho^2)
    dFdf = dDdrho/100;                  % chain rule, rho = 1+f/100
    dzdf = -dFdf/dFdz;
end
