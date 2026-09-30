%% Part 5: Root-locus departure sensitivity - approximate K_max(f)
%
% Objective: derive a first-order (small-K) analytical explanation for
% the results of Part 4 using classical root-locus departure theory,
% and check it against the numerically exact stability boundary.
%
% For 1 + K*G_OL(z) = 0, i.e. D(z,f) + K*N(z) = 0, the velocity at
% which a closed-loop root leaves the open-loop pole p as K increases
% from 0 is the standard result
%       dz/dK |_{K=0} = -N(p) / D'(p,f)
% Evaluating this at the complex open-loop pole p(f) = rho*exp(j*2pi/3)
% tells us, to first order, whether positive gain pushes that pole
% into or out of the unit circle, and (by requiring |p+K*dz/dK| = 1)
% gives an approximate closed-form stability boundary K_max(f).

num_GOL = [1 1.7 -1];

fprintf('%6s  %-18s  %14s  %s\n','f(%)','p(f)','d|z|^2/dK','direction')
for f = [-30 -20 -10 -5 -2.5 -1 0 1 5 10 20 30]
    rho = 1 + f/100;
    p = rho*exp(1j*2*pi/3);
    den_f = conv([1 -1],[1 rho rho^2]);
    Dp = polyval(polyder(den_f), p);
    Np = polyval(num_GOL, p);
    dzdK  = -Np/Dp;                     % departure velocity
    slope = 2*real(conj(p)*dzdK);       % d|z|^2/dK at K=0
    if slope < 0, dirn = 'INWARD (stabilising)';
    else,          dirn = 'OUTWARD (destabilising)'; end
    fprintf('%6.1f  %8.4f%+8.4fi  %14.5f  %s\n', f, real(p), imag(p), slope, dirn)
end
disp('=> the complex pole departs OUTWARD for every f tested: any K>0 makes')
disp('   this mode worse from the outset, which is why f>=0 is never stabilisable.')

% ---- First-order approximate stability boundary K_max(f), for f<0 ----
fvec = linspace(-30,-0.1,60);
K_boundary_linearised = nan(size(fvec));
for i = 1:numel(fvec)
    f = fvec(i);
    rho = 1+f/100;
    p = rho*exp(1j*2*pi/3);
    den_f = conv([1 -1],[1 rho rho^2]);
    Dp = polyval(polyder(den_f), p);
    Np = polyval(num_GOL, p);
    dzdK = -Np/Dp;
    A = abs(dzdK)^2;
    B = 2*real(conj(p)*dzdK);
    C = abs(p)^2 - 1;
    r = roots([A B C]);
    r = r(abs(imag(r))<1e-9 & real(r)>0);
    if ~isempty(r)
        K_boundary_linearised(i) = min(real(r));
    end
end
disp('K_boundary_linearised(f) computed (compare against the numerically exact')
disp('K_max(f) from Part 4''s stability scan - see Fig. 7 in the report).')
