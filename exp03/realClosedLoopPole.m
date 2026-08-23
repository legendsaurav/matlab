function s3 = realClosedLoopPole(a, K)
%REALCLOSEDLOOPPOLE  The (most nearly) real closed-loop pole.
%   s3 = REALCLOSEDLOOPPOLE(a, K) returns the root of
%   s^3+a*s^2+K*s+K=0 with the smallest |Im|, i.e. the pole that tracks
%   toward the open-loop zero at s=-1 as K -> Inf (see report Figure 12).
    r = roots(charEqCoeffs(a, K));
    [~, idx] = min(abs(imag(r)));
    s3 = real(r(idx));
end
