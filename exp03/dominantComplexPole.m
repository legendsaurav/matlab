function s = dominantComplexPole(a, K)
%DOMINANTCOMPLEXPOLE  Dominant (largest |Im|) closed-loop pole.
%   s = DOMINANTCOMPLEXPOLE(a, K) solves s^3+a*s^2+K*s+K=0 and returns the
%   root with the largest imaginary magnitude. This is a robust
%   generalisation of the original lab code's fixed threshold
%   "poles(imag(poles) > 0.1)": that threshold only works when K is large
%   (imag part ~ sqrt(K) >> 0.1); selecting by max(|Im|) instead degrades
%   gracefully for the small-K sweeps used in Part D.
    r = roots(charEqCoeffs(a, K));
    [~, idx] = max(abs(imag(r)));
    s = r(idx);
end
