function Kp_stab = marginal_stability_gain(b0,p)
%MARGINAL_STABILITY_GAIN  Exact gain at which a closed-loop pole first
%reaches the unit circle (the true stability limit).
%
%   Kp_stab = marginal_stability_gain(b0,p)
%
%   On the complex branch, |z|^2 equals the product of the roots of
%   z^2 + B z + C = 0, which is the constant term C = p + Kp*b0. The loop
%   goes unstable when that product reaches 1:
%
%       Kp_stab = (1 - p) / b0
%
%   See EE208 Experiment 7 report, Sec. 4.7 / 6.1.

    Kp_stab = (1 - p)/b0;
end
