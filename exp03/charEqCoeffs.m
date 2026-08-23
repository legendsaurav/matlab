function p = charEqCoeffs(a, K)
%CHAREQCOEFFS  Coefficients of s^3 + a*s^2 + K*s + K = 0
%   p = CHAREQCOEFFS(a, K) returns the coefficient row vector [1 a K K]
%   suitable for MATLAB's roots().
    p = [1, a, K, K];
end
