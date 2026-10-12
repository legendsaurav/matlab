function [Kp0,c,R2] = fit_inverse_law(Tvec,Kpvec)
%FIT_INVERSE_LAW  Fit Kp(T) ~= Kp0 / (1 + c*T) by nonlinear least squares.
%
%   [Kp0,c,R2] = fit_inverse_law(Tvec,Kpvec)
%
%   Uses FMINSEARCH (base MATLAB) instead of lsqcurvefit, so there is no
%   Optimization Toolbox dependency.
%   See EE208 Experiment 7 report, Sec. 4.8 / 6.3.

    Tvec = Tvec(:);  Kpvec = Kpvec(:);
    model = @(par,T) par(1)./(1+par(2).*T);
    sse   = @(par) sum((Kpvec - model(par,Tvec)).^2);

    par0 = [Kpvec(1), 1];
    opts = optimset('Display','off','TolX',1e-12,'TolFun',1e-14, ...
                    'MaxFunEvals',5000,'MaxIter',5000);
    par  = fminsearch(sse, par0, opts);

    Kp0 = par(1);
    c   = par(2);
    resid = Kpvec - model(par,Tvec);
    R2 = 1 - sum(resid.^2)/sum((Kpvec-mean(Kpvec)).^2);
end
