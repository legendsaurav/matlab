function part2_rootlocus_stability()
%% Part 2: Root locus and stability sweep of the closed loop vs gain K
%
% Objective: trace how the three closed-loop poles move in the
% z-plane as the proportional gain K varies (the discrete root
% locus of 1 + K*G_OL(z) = 0), and determine numerically whether any
% real value of K places all three closed-loop poles strictly inside
% the unit circle (i.e. whether proportional control alone can
% stabilise this plant).

num_GOL = [1 1.7 -1];
den_GOL = [1 0 0 -1];

% Closed-loop characteristic polynomial for a given gain K
closedLoopPoly = @(K) den_GOL + K*[0 num_GOL];

% ---- Root-locus branches (K>0 and K<0), tracked for continuity ----
Kpos = linspace(0, 60, 3000);
Kneg = linspace(0, -60, 3000);
rlocusPos = trackBranches(Kpos, closedLoopPoly);
rlocusNeg = trackBranches(Kneg, closedLoopPoly);
% (equivalently, rlocus(tf(num_GOL,den_GOL,Ts)) from the Control
%  System Toolbox produces the same branches for K>0)

% ---- Stability sweep: largest closed-loop pole magnitude vs K ----
Ksweep = linspace(-15, 15, 601);
maxAbsRoot = zeros(size(Ksweep));
for i = 1:numel(Ksweep)
    maxAbsRoot(i) = max(abs(roots(closedLoopPoly(Ksweep(i)))));
end
isStable = maxAbsRoot < 1;
fprintf('Any K in [-15,15] with all closed-loop poles strictly inside unit circle? %d\n', any(isStable))
[minMag, iMin] = min(maxAbsRoot);
fprintf('min(max|root(K)|) = %.6f, attained at K = %.4f\n', minMag, Ksweep(iMin))

% Wider confirmation (K up to +-1000, log-spaced away from 0)
Kwide = [-fliplr(logspace(-2,6,4000)), logspace(-2,6,4000)];
maxAbsRootWide = zeros(size(Kwide));
for i = 1:numel(Kwide)
    maxAbsRootWide(i) = max(abs(roots(closedLoopPoly(Kwide(i)))));
end
fprintf('Global min(max|root|) over K in +-[1e-2,1e6]: %.6f\n', min(maxAbsRootWide))
disp('=> No real proportional gain K stabilises this plant (confirms Sec. 6 discussion).')
end

%% ---- local functions ----
function rootsSorted = trackBranches(Ks, polyFcn)
% Track the 3 characteristic-equation roots continuously as K is
% swept, by greedily matching each new root to its nearest neighbour
% from the previous step (standard trick for drawing coherent root
% locus branches from a bare root-finder).
    n = numel(Ks);
    rootsSorted = zeros(n,3);
    prev = [];
    for i = 1:n
        r = roots(polyFcn(Ks(i)));
        if ~isempty(prev)
            used = false(1,3);
            order = zeros(1,3);
            for a = 1:3
                d = abs(prev(a) - r);
                d(used) = inf;
                [~,j] = min(d);
                used(j) = true;
                order(a) = j;
            end
            r = r(order);
        end
        rootsSorted(i,:) = r;
        prev = r;
    end
end
