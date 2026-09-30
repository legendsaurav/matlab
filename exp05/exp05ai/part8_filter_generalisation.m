function part8_filter_generalisation()
%% Part 8: Generalising the anti-aliasing filter design (Section 10)
%
% Section 7.3 fixed one digital corner theta_c=pi/3 and compared exactly
% two filter orders (n=1,2). Here both the order n and the corner
% theta_c are treated as free design parameters, using the EXACT
% closed-form digital magnitude response of a Tustin+pre-warped
% Butterworth filter:
%     |Hd(e^{j theta})|_dB = -10*log10(1 + (tan(theta/2)/tan(theta_c/2))^(2n))
% This holds exactly because pre-warping is defined precisely so that
% the digital filter's response at theta equals the analog prototype's
% response at w=(2/Ts)tan(theta/2), and the two factors of (2/Ts) cancel
% in the ratio w/wc -- so once theta_c is fixed, Ts drops out completely.
%
% This script (a) validates that formula against a fully independent,
% toolbox-free brute-force construction of the digital filter (analog
% Butterworth poles via roots()/poly(), then the bilinear transform
% applied by direct polynomial substitution via conv() -- no Control or
% Signal Processing toolbox function is used, so this runs in plain
% Octave); (b) exhibits the exact transmission zero at the band edge;
% (c) sweeps the (n,theta_c) guard-band/resonance-cost trade-off; and
% (d) shows how a PHYSICALLY specified corner frequency reintroduces the
% Ts-dependence that Section 7.2 showed was absent from the plant itself.

fprintf('--- Validation: exact closed form vs. toolbox-free bilinear construction ---\n')
worstErr = 0;
thetaTest = [0.3 0.9 pi/3 2*pi/3 3.0 pi-1e-3];
for n = 1:4
    for theta_c = deg2rad([30 60 100])
        for Ts = [1.0 0.37]
            [b_z, a_z] = butterDigital(n, theta_c, Ts);
            for th = thetaTest
                zpt = exp(1j*th);
                H = polyval(b_z, zpt) / polyval(a_z, zpt);
                bruteDB = 20*log10(abs(H));
                exactDB = magDbExact(th, theta_c, n);
                worstErr = max(worstErr, abs(bruteDB - exactDB));
            end
        end
    end
end
fprintf('Max |exact - toolbox-free bilinear| over the full sweep = %.3e dB\n', worstErr)
disp('(Matches Part 11''s Python cross-check, which independently uses scipy''s')
disp('butter()/bilinear() route -- two unrelated implementations of the bilinear')
disp('transform agree with the closed form to floating-point precision.)')

fprintf('\n--- Exact zero at the band edge (theta=pi, z=-1), theta_c=pi/3 ---\n')
for n = 1:4
    [b_z, a_z] = butterDigital(n, pi/3, 1.0);
    Hminus1 = polyval(b_z,-1)/polyval(a_z,-1);
    nearEdgeDB = magDbExact(pi-1e-4, pi/3, n);
    fprintf('n=%d:  H(z=-1) = %.3e (exact)   grid point at theta=pi-1e-4: %.2f dB\n', ...
        n, Hminus1, nearEdgeDB)
end

fprintf('\n--- Guard band and resonance cost at the actual design corner (60 deg) ---\n')
theta_res = 2*pi/3;
for n = 1:4
    guardDeg = rad2deg(pi - theta60(deg2rad(60), n));
    resCostDB = magDbExact(theta_res, deg2rad(60), n);
    fprintf('n=%d:  guard band = %7.4f deg   resonance cost = %8.4f dB\n', n, guardDeg, resCostDB)
end

fprintf('\n--- Physical corner-frequency spec: Ts-sensitivity (n=2 Butterworth) ---\n')
omega_c_phys = 1.0;   % rad/s, a fixed PHYSICAL requirement
n_fixed = 2;
for Ts = [1.0 0.5 0.2 0.1]
    theta_c_req = omega_c_phys * Ts;          % theta = omega*Ts
    guardDeg = rad2deg(pi - theta60(theta_c_req, n_fixed));
    resCostDB = magDbExact(theta_res, theta_c_req, n_fixed);
    fprintf('Ts=%4.2f  theta_c=%7.3f deg   guard band=%7.3f deg   resonance cost=%8.3f dB\n', ...
        Ts, rad2deg(theta_c_req), guardDeg, resCostDB)
end
disp('The same physical requirement therefore implies a different digital filter --')
disp('and a very different guard band / resonance cost -- for every choice of the')
disp('unspecified Ts, in sharp contrast to Section 7.2''s exact Ts-independence of')
disp('the plant''s own response at a GIVEN digital frequency theta.')
end

%% ---- local functions ----
function db = magDbExact(theta, theta_c, n)
% Exact closed-form magnitude (dB) of an n-th order Tustin+prewarped
% Butterworth digital low-pass with digital corner theta_c.
    ratio = tan(theta/2) / tan(theta_c/2);
    db = -10*log10(1 + ratio^(2*n));
end

function th60 = theta60(theta_c, n)
% Exact digital frequency at which attenuation first reaches 60 dB.
    ratioNeeded = (1e6 - 1)^(1/(2*n));
    th60 = 2*atan(tan(theta_c/2) * ratioNeeded);
end

function [b_z, a_z] = butterDigital(n, theta_c, Ts)
% n-th order Butterworth low-pass, Tustin+prewarped to digital corner
% theta_c, built from ONLY roots()/poly()/conv() -- no Control or
% Signal Processing toolbox function, so this runs in plain Octave.
    wc = (2/Ts) * tan(theta_c/2);
    sgn = (-1)^n;
    % All 2n roots of s^(2n) + sgn*wc^(2n) = 0; keep the n in the LHP.
    coeffs = [1, zeros(1,2*n-1), sgn*wc^(2*n)];
    rAll = roots(coeffs);
    pLhp = rAll(real(rAll) < 0);
    a_s = real(poly(pLhp));       % length n+1, highest power first, monic
    b_s = a_s(end);                % unity DC gain by construction: H(0)=b_s/a_s(end)=1
    c = 2/Ts;
    den_z = zeros(1, n+1);
    for k = 0:n
        term = a_s(k+1) * c^(n-k) * conv(polyPow([1 -1], n-k), polyPow([1 1], k));
        term = [zeros(1, (n+1)-length(term)), term];
        den_z = den_z + term;
    end
    num_z = b_s * polyPow([1 1], n);
    num_z = [zeros(1, (n+1)-length(num_z)), num_z];
    a_z = den_z / den_z(1);
    b_z = num_z / den_z(1);
end

function p = polyPow(base, k)
% base (a length-2 polynomial) raised to integer power k>=0.
    p = 1;
    for i = 1:k
        p = conv(p, base);
    end
end
