function part4_overshoot_sensitivity()
%% Part 4: Overshoot sensitivity to gain K and pole-radius perturbation f
%
% Objective: since Part 2 shows the nominal plant cannot be stabilised
% by any proportional gain, ask a natural follow-up question: if the
% open-loop resonant pole pair were only known to +-30% in radius (a
% plausible identification/modelling uncertainty, parameterised by
% f, in percent), how would the closed-loop step response's Percent
% Overshoot (PO) depend on K and f jointly, and can PO be predicted
% by a simple approximate formula?
%
% Modelling choice: only the complex pole pair rho*exp(+-j*2*pi/3) is
% perturbed (rho = 1+f/100); the integrator pole at z=1 and the zeros
% of G_OL(z) are structurally fixed by the given characteristic
% polynomial, so they are left unperturbed. f = 0 reproduces the
% original plant exactly.

num_GOL = [1 1.7 -1];

perturbedDen = @(f) conv([1 -1], [1 (1+f/100) (1+f/100)^2]);
closedLoopPolyF = @(K,den_f) den_f + K*[0 num_GOL];

% ---- Worked examples at K = 0.05 for f = -20, 0, +20 %  ----
for f = [-20 0 20]
    K = 0.05;
    r = roots(closedLoopPolyF(K, perturbedDen(f)));
    [~, idx] = max(abs(r));
    zd = r(idx);
    y  = stepResponseClosedLoop(K, f, 150);
    PO_exact  = percentOvershoot(y, 1.0, 100);
    PO_approx = poApproxFromPole(zd);
    fprintf(['K=%.2f  f=%+4d%%  dominant pole=%6.3f%+6.3fi  stable=%d  ', ...
             'PO_exact=%6.1f%%  PO_approx=%6.1f%%\n'], ...
        K, f, real(zd), imag(zd), max(abs(r))<1, PO_exact, PO_approx)
end

% The full sensitivity map PO(K,f) used in the report (Fig. 8-9) is
% obtained by evaluating stepResponseClosedLoop/percentOvershoot and
% poApproxFromPole over a grid of (K,f) pairs; see
% run_full_analysis.m for the grid-sweep and plotting code.
end

%% ---- local functions ----
function y = stepResponseClosedLoop(K, f_percent, n)
% Unit-step response of the unity-feedback closed loop, obtained by
% direct recursion on the difference equation implied by
%   T(z) = K*N(z) / ( D(z,f) + K*N(z) )
% A short horizon and a per-sample amplitude clip keep the result
% finite and smooth even when (K,f) lies in the unstable region.
    num_GOL = [1 1.7 -1];
    rho = 1 + f_percent/100;
    den_f = conv([1 -1],[1 rho rho^2]);
    a = den_f + K*[0 num_GOL];
    b = K*num_GOL;
    a = a/a(1);
    bpad = [0 b];
    clipVal = 1e6;
    y = zeros(1,n); u = ones(1,n);
    for k = 1:n
        acc = 0;
        for j = 2:4
            if k-j+1 >= 1
                acc = acc + bpad(j)*u(k-j+1) - a(j)*y(k-j+1);
            end
        end
        yk = acc/a(1);
        if ~isfinite(yk); yk = clipVal; end
        y(k) = max(min(yk,clipVal),-clipVal);
    end
end

function po = percentOvershoot(y, finalValue, cap)
% Finite-horizon, saturating percent-overshoot metric:
%   PO = 100*(max(y)-finalValue)/finalValue, clipped to [0,cap]
    y(~isfinite(y)) = cap*10;
    po = 100*(max(y)-finalValue)/finalValue;
    po = min(max(po,0),cap);
end

function po = poApproxFromPole(zd)
% Dominant-pole approximate percent overshoot from a single discrete
% closed-loop pole z_d = r*exp(j*phi):
%   zeta = -ln(r) / sqrt( ln(r)^2 + phi^2 )      (equivalent damping ratio)
%   PO   = 100 * exp( pi*ln(r)/phi )             (classical 2nd-order formula)
% Both expressions are independent of the sampling time Ts.
    r = abs(zd); phi = abs(angle(zd));
    if phi < 1e-9 || r <= 0
        po = 0; return
    end
    lr = log(r);
    zeta = -lr/sqrt(lr^2+phi^2); %#ok<NASGU>  (kept for reporting/inspection)
    po = 100*exp(pi*lr/phi);
    po = min(max(po,0),200);
end
