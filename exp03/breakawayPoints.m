function [D, s_bi, K_bi, s_ba, K_ba] = breakawayPoints(a)
%BREAKAWAYPOINTS  Breakaway/break-in analysis for a fixed a.
%   [D, s_bi, K_bi, s_ba, K_ba] = BREAKAWAYPOINTS(a)
%   D      : discriminant (a+7.5)^2-40a = (a-2.5)(a-22.5)
%   s_bi   : break-in point (real s, the LOWER-gain event)   [NaN if D<0]
%   K_bi   : gain at break-in
%   s_ba   : breakaway point (real s, the HIGHER-gain event) [NaN if D<0]
%   K_ba   : gain at breakaway
%
%   Derived from dK/ds = 0 on K(s) = -s^2(s+a)/(s+2.5), which factors as
%   s*[2s^2+(a+7.5)s+5a] = 0. The quadratic factor's two roots are the
%   candidates handled here (s=0, the trivial origin breakaway, is not
%   returned since it exists for every a).
    D = (a + 7.5)^2 - 40*a;
    if D < 0
        s_bi = NaN; K_bi = NaN; s_ba = NaN; K_ba = NaN;
        return
    end
    s1 = (-(a+7.5) + sqrt(D)) / 4;
    s2 = (-(a+7.5) - sqrt(D)) / 4;
    K1 = -s1^2*(s1+a) / (s1+2.5);
    K2 = -s2^2*(s2+a) / (s2+2.5);
    if K1 <= K2
        s_bi = s1; K_bi = K1; s_ba = s2; K_ba = K2;
    else
        s_bi = s2; K_bi = K2; s_ba = s1; K_ba = K1;
    end
end
