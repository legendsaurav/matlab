function Kp_target = zeta_crossing_gain(b1,b0,p,target_zeta,Kp_lo,Kp_hi)
%ZETA_CROSSING_GAIN  Gain at which the complex-pole branch crosses a
%constant-damping-ratio contour (e.g. zeta = 0.6901 for 5% overshoot).
%
%   Kp_target = zeta_crossing_gain(b1,b0,p,target_zeta,Kp_lo,Kp_hi)
%
%   Uses the T-independent EFFECTIVE_ZETA, so the same routine works at
%   every sampling time. [Kp_lo,Kp_hi] must bracket the root: use the two
%   outputs of BREAKAWAY_GAIN, nudged slightly inward
%   (Kp_min*1.0005, Kp_breakin*0.999) to avoid the repeated root at each
%   end.
%
%   FZERO is base MATLAB (no toolbox needed).
%   See EE208 Experiment 7 report, Sec. 4.5.

    Kp_target = fzero(@(Kp) zeta_error(Kp,b1,b0,p,target_zeta), [Kp_lo Kp_hi]);
end
