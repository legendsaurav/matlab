function [Kp_min,Kp_breakin] = breakaway_gain(b1,b0,p)
%BREAKAWAY_GAIN  Exact zeta = 1 boundary of the digital root locus.
%
%   [Kp_min,Kp_breakin] = breakaway_gain(b1,b0,p)
%
%   Closed-loop characteristic equation:
%       (z-1)(z-p) + Kp*(b1*z+b0) = 0
%       z^2 + [Kp*b1-(1+p)] z + [p+Kp*b0] = 0
%
%   The locus leaves the real axis (zeta = 1) where the two roots first
%   coincide, i.e. where the discriminant vanishes. That condition is a
%   quadratic in Kp:
%
%       b1^2*Kp^2 - [2*b1*(1+p) + 4*b0]*Kp + (1-p)^2 = 0
%
%   Its two roots are the breakaway gain (poles turn complex) and the
%   break-in gain (poles rejoin the real axis). Kp_breakin is returned
%   as the upper bracket for ZETA_CROSSING_GAIN.
%
%   See EE208 Experiment 7 report, Sec. 4.4.

    A = b1^2;
    B = -(2*b1*(1+p) + 4*b0);
    C = (1-p)^2;
    r = sort(real(roots([A B C])));
    Kp_min     = r(1);
    Kp_breakin = r(2);
end
