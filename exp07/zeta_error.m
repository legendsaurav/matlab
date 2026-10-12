function e = zeta_error(Kp,b1,b0,p,target_zeta)
%ZETA_ERROR  Dominant-pole damping ratio minus a target, as a function of
%gain. This is the function FZERO drives to zero in ZETA_CROSSING_GAIN
%and in the lead-compensator search (Sec. 13), where p is the
%compensator's free pole p_lead instead of the plant pole.
%
%   e = zeta_error(Kp,b1,b0,p,target_zeta)

    B  = Kp*b1 - (1+p);
    C  = p + Kp*b0;
    r  = roots([1 B C]);
    zc = r(imag(r) > 1e-12);       % upper complex-conjugate pole
    if isempty(zc)
        % Below breakaway the poles are still real: no complex pole to
        % measure. Treat as "more damped than any target".
        e = 1 - target_zeta;
        return
    end
    e = effective_zeta(zc(1)) - target_zeta;
end
