function [b1,b0,p] = zoh_coeffs(T,a,K)
%ZOH_COEFFS  ZOH-equivalent numerator coefficients of K/[s(s+a)] at
%sampling time T, in closed form.
%
%   [b1,b0,p] = zoh_coeffs(T,a,K) returns the ZOH discretisation
%       G(z) = (b1*z + b0) / [(z-1)(z-p)]
%   of the continuous plant K/[s(s+a)]. Equivalent to
%       Gz = c2d(tf(K,[1 a 0]), T, 'zoh');
%   but in closed form, so the root-finding functions in this package
%   can use b1, b0, p directly. T may be a scalar or a vector.
%
%   See EE208 Experiment 7 report, Sec. 3.3 / 4.3.

    p  = exp(-a*T);
    b1 = (K/a^2)*(a*T - 1 + p);
    b0 = (K/a^2)*(1 - p - a*p.*T);
end
