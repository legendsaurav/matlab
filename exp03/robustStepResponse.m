function [t, y] = robustStepResponse(a, K, TendOverride)
%ROBUSTSTEPRESPONSE  True step response of T(s)=K(s+1)/(s^3+as^2+Ks+K).
%   [t, y] = ROBUSTSTEPRESPONSE(a, K) simulates the unit-step response
%   using an EXPLICIT, fine time grid built from the raw closed-loop
%   pole magnitudes -- NOT a zeta/Tpeak trig formula -- so it degrades
%   gracefully whether the dominant pole is complex OR real.
%
%   *** REVISION NOTE ***
%   An earlier version of this function estimated the required time
%   window from Tpeak = pi/(wn*sqrt(1-zeta^2)). That blows up as
%   zeta -> 1, which is exactly what happens for a > 9 at gains INSIDE
%   the breakaway/break-in loop of Sec. 5.3 (all three closed-loop poles
%   real, so the "dominant complex pole" selector picks a real pole with
%   near-zero imaginary part, giving zeta ~ 1). That produced Tend in
%   the tens of thousands of seconds and a 4000-point grid spread across
%   it -- so coarse that the entire early transient, including the true
%   peak, was missed, silently returning ~0% overshoot instead of the
%   correct value. Example: a=20, K=90 returned 0.00% instead of the
%   correct 15.23%.
%
%   This version instead uses:
%     decay_min = min(|Re(poles)|)   -- slowest pole, sets settling time
%     w_max     = max(|poles|)        -- fastest pole/frequency, sets
%                                         the resolution needed
%   Both are well-defined and well-behaved for real OR complex poles,
%   so there is no failure mode analogous to the one above. Verified
%   against fine-grid (300k-3M point) reference simulations across
%   lightly-damped, heavily-damped, and all-real-pole regimes -- see
%   Sec. 10.4 of the report for the a=20 case that exposed the bug.
%
%   [t, y] = ROBUSTSTEPRESPONSE(a, K, TendOverride) uses a caller-supplied
%   final time instead of the automatically-estimated one (e.g. to match
%   a specific plot window).
    r = roots(charEqCoeffs(a, K));
    decay_min = min(abs(real(r)));    % slowest pole -> settling time
    w_max = max(abs(r));               % fastest pole/freq -> resolution

    Ts_est = 5 / max(decay_min, 1e-6);

    if nargin >= 3 && ~isempty(TendOverride)
        Tend = TendOverride;
    else
        Tend = max([1.3*Ts_est, 10/w_max]);
    end

    dt_needed = (1/w_max) / 50;
    N = max(4000, ceil(Tend / dt_needed));
    N = min(N, 2e5);   % sane upper cap on point count

    tgrid = linspace(0, Tend, N);
    sys = tf([K K], charEqCoeffs(a, K));
    [y, t] = step(sys, tgrid);
end
