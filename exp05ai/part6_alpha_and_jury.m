function part6_alpha_and_jury()
%% Part 6: Generalised structural coefficient alpha - exact stability region
%
% Objective: Part 1-5 fix the numerator coefficient in G_OL(z) at its
% given value, 1.7. Here that coefficient is generalised to a free
% parameter alpha (alpha=1.7 recovers the original problem exactly),
% and Jury's stability criterion (the discrete-time analogue of
% Routh-Hurwitz) is applied SYMBOLICALLY to the closed-loop polynomial
%     z^3 + K*z^2 + alpha*K*z - (K+1) = 0
% to determine, in closed form, every (alpha,K) pair for which the
% closed loop is stable. This both (i) turns the earlier ad-hoc
% numerical stability sweep into a proved, closed-form result, and
% (ii) places the actual plant (alpha=1.7) precisely relative to the
% region of plants that COULD have been stabilised.
%
% Jury conditions for a monic cubic z^3+a2 z^2+a1 z+a0 (all three
% roots inside the unit circle) are P(1)>0, -P(-1)>0, |a0|<1 and
% |b0|>|b2| with b0=a0^2-1, b2=a0*a2-a1. Substituting a2=K, a1=alpha*K,
% a0=-(K+1) and simplifying (checked independently with symbolic
% algebra) gives the four scalar conditions implemented below.

fprintf('--- Jury conditions at the actual plant, alpha=1.7 ---\n')
for K = [-1.5 -0.5 0.1 0.5 1.5]
    [ok,j1,j2,j3,j4] = juryStable(1.7, K);
    fprintf('K=%6.2f   J1=%d J2=%d J3=%d J4=%d   all-stable=%d\n', K, j1,j2,j3,j4, ok)
end
disp('J1 (alpha*K>0) and J3 (|K+1|<1, i.e. -2<K<0) can never both hold when')
disp('alpha=1.7>0, so NO real K stabilises the actual plant -- confirming the')
disp('numerical stability sweep of Part 2 with a closed-form proof.')

% ---- Exact closed-form region: alpha in (-3,0),  K in (-(3+alpha)/2, 0) ----
fprintf('\n--- Exact stability window K in (K_lower(alpha), 0) ---\n')
for alpha = [-2.9 -2.5 -2.0 -1.0 -0.5 -0.1]
    Klo = alphaKLowerBound(alpha);
    fprintf('alpha=%6.2f   K_lower = -(3+alpha)/2 = %7.4f\n', alpha, Klo)
end

% ---- Optimal design point: triple coalesced closed-loop pole ----
% At the optimum all three closed-loop poles coincide at a single real
% value z=r*, which forces (matching coefficients with (z-r)^3):
%   K* = -3r*,  alpha* = -r*,  K* = r*^3 - 1
% Eliminating K* gives the depressed cubic  r^3 + 3r - 1 = 0, solved
% here in closed form via Cardano's formula (the golden ratio phi
% appears inside the cube roots -- a curiosity of this particular
% cubic, not a general control-theory fact).
phi = (1+sqrt(5))/2;
r_star = phi^(1/3) - phi^(-1/3);
alpha_star = -r_star;
K_star = r_star^3 - 1;
p_star = [1 K_star alpha_star*K_star -(K_star+1)];
roots_star = roots(p_star);
fprintf('\n--- Optimal (alpha*,K*): minimises the dominant closed-loop pole magnitude ---\n')
fprintf('r*=%.10f  alpha*=%.10f  K*=%.10f\n', r_star, alpha_star, K_star)
fprintf('closed-loop roots at optimum: '); disp(roots_star.')
fprintf('max|root| = %.10f (equals r* to machine precision)\n', max(abs(roots_star)))

% ---- Exact K_max(f) at the nominal plant (alpha=1.7), from the fully
% general 3-parameter (K,rho,alpha) Jury condition (J4) specialised to
% alpha=1.7 and solved for K(rho); rho=1+f/100. Compare with the
% first-order (linearised) approximation of Part 5. ----
fprintf('\n--- Exact vs. linearised K_max(f), alpha=1.7 ---\n')
for f = [-30 -20 -10 -5 -2.5]
    Kex = KmaxExactJury(f);
    fprintf('f=%6.1f%%   K_max (exact, Jury) = %.5f\n', f, Kex)
end
disp('See table7_kmax_exact_vs_linearised.csv / Fig. 15 for the full curve and')
disp('its comparison against Part 5''s first-order departure-velocity estimate.')
end

%% ---- local functions ----
function [allStable,j1,j2,j3,j4] = juryStable(alpha, K)
% Jury stability test for z^3 + K z^2 + alpha*K z - (K+1) = 0.
    a2 = K; a1 = alpha*K; a0 = -(K+1);
    P1  = 1 + a2 + a1 + a0;      % = alpha*K after simplification
    nPm1 = 1 - a2 + a1 - a0;     % = -P(-1)
    b0 = a0^2 - 1;
    b2 = a0*a2 - a1;
    j1 = P1 > 0;
    j2 = nPm1 > 0;
    j3 = abs(a0) < 1;
    j4 = (b0-b2)*(b0+b2) > 0;
    allStable = j1 && j2 && j3 && j4;
end

function Klo = alphaKLowerBound(alpha)
% Exact lower edge of the (alpha,K) stability region, valid for
% -3 < alpha < 0 (upper edge is K=0 for all such alpha).
    Klo = -(3+alpha)/2;
end

function K = KmaxExactJury(f_percent)
% Exact stabilising-gain boundary at alpha=1.7, obtained by solving
% the general 3-parameter Jury condition (J4) for K at fixed
% rho=1+f/100 (symbolic derivation carried out with computer algebra;
% reproduced here as a closed-form expression).
    rho = 1 + f_percent/100;
    disc = 100*rho^4 - 200*rho^3 + 520*rho^2 + 940*rho + 849;
    K = -0.75*rho^2 - 0.25*rho + sqrt(disc)/40 - 7/40;
end
