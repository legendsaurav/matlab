function [D, s_bi, K_bi, s_ba, K_ba] = breakawayPoints(a)
%BREAKAWAYPOINTS  Breakaway/break-in analysis for a fixed a.
%   [D, s_bi, K_bi, s_ba, K_ba] = BREAKAWAYPOINTS(a)
%   D      : discriminant (a+3)^2-16a = (a-1)(a-9)
%   s_bi   : break-in point (real s, the LOWER-gain event)   [NaN if D<0]
%   K_bi   : gain at break-in
%   s_ba   : breakaway point (real s, the HIGHER-gain event) [NaN if D<0]
%   K_ba   : gain at breakaway
%
%   Derived from dK/ds = 0 on K(s) = -s^2(s+a)/(s+1), which factors as
%   s*[2s^2+(a+3)s+2a] = 0. The quadratic factor's two roots are the
%   candidates handled here (s=0, the trivial origin breakaway, is not
%   returned since it exists for every a).
    D = (a + 3)^2 - 16*a;
    if D < 0
        s_bi = NaN; K_bi = NaN; s_ba = NaN; K_ba = NaN;
        return
    end
    s1 = (-(a+3) + sqrt(D)) / 4;
    s2 = (-(a+3) - sqrt(D)) / 4;
    K1 = -s1^2*(s1+a) / (s1+1);
    K2 = -s2^2*(s2+a) / (s2+1);
    if K1 <= K2
        s_bi = s1; K_bi = K1; s_ba = s2; K_ba = K2;
    else
        s_bi = s2; K_bi = K2; s_ba = s1; K_ba = K1;
    end
end
