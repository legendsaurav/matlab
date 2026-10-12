function [dzdKp,zc] = pole_sensitivity(Kp,b1,b0,p)
%POLE_SENSITIVITY  dz/dKp of the dominant closed-loop pole, by implicit
%differentiation of the characteristic equation.
%
%   [dzdKp,zc] = pole_sensitivity(Kp,b1,b0,p)
%
%   F(z,Kp) = z^2 + [Kp*b1-(1+p)] z + [p+Kp*b0] = 0  gives
%
%       dz/dKp = -(b1*z + b0) / (2z + Kp*b1 - (1+p))
%
%   The denominator is the breakaway discriminant condition, so the
%   sensitivity diverges at breakaway and shrinks farther along the
%   locus. Also returns the pole zc used (upper half plane if complex).
%
%   See EE208 Experiment 7 report, Sec. 12.1.

    B = Kp*b1 - (1+p);
    C = p + Kp*b0;
    r = roots([1 B C]);
    if any(abs(imag(r)) > 1e-9)
        zc = r(imag(r) > 0);
        zc = zc(1);
    else
        zc = max(r);               % real pole (Kp below breakaway)
    end
    dzdKp = -(b1*zc + b0)/(2*zc + B);
end
