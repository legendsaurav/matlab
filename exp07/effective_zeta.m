function zeta = effective_zeta(z)
%EFFECTIVE_ZETA  Damping ratio of a z-domain pole, independent of T.
%
%   zeta = effective_zeta(z) maps a (possibly complex) closed-loop pole z
%   back to the damping ratio of the continuous pole s = ln(z)/T. The
%   result does not depend on T: changing T only rescales |s|, never its
%   angle, and zeta depends on the angle alone.
%
%       zeta_eff(z) = -ln|z| / sqrt( (ln|z|)^2 + (arg z)^2 )
%
%   Elementwise: works on a scalar, vector or matrix of poles (uses ./,
%   because / between same-size vectors is matrix division in MATLAB).
%
%   See EE208 Experiment 7 report, Sec. 3.3.

    lz   = log(z);                 % complex log: ln|z| + j*arg(z)
    zeta = -real(lz)./abs(lz);
end
