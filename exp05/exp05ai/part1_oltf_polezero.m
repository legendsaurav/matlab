%% Part 1: Recover the open-loop transfer function; analyse poles & zeros
%
% Objective: given only the closed-loop characteristic polynomial
%   z^3 + K z^2 + 1.7 K z - (K+1) = 0
% recover the open-loop transfer function G_OL(z) of the unity-
% negative-feedback digital plant, then locate and interpret its
% poles and zeros in the z-plane.
%
% Method: split the polynomial into the part independent of K and the
% part multiplying K, then match against the standard unity-feedback
% closed-loop characteristic equation 1 + K*G_OL(z) = 0:
%
%   (z^3 - 1)  +  K (z^2 + 1.7z - 1)  =  0
%   1 + K * (z^2+1.7z-1)/(z^3-1)      =  0     =>  G_OL(z) = (z^2+1.7z-1)/(z^3-1)

num_GOL = [1 1.7 -1];      % numerator:   z^2 + 1.7z - 1  (zeros)
den_GOL = [1 0 0 -1];      % denominator: z^3 - 1         (poles)

zz = roots(num_GOL);       % open-loop zeros
pp = roots(den_GOL);       % open-loop poles

disp('--- Zeros of G_OL(z) ---'); disp(zz)
disp('--- Poles of G_OL(z) ---'); disp(pp)
disp('--- Pole magnitudes ---');  disp(abs(pp))
disp('--- Pole angles (deg) ---'); disp(angle(pp)*180/pi)

% Verify the recovered G_OL(z) actually reproduces the given
% closed-loop polynomial for arbitrary K (algebraic sanity check)
for K = [-3.3 0 1 4.2 7.75]
    lhs = den_GOL + K*[0 1 1.7 -1];      % den_GOL + K*num_GOL, aligned
    rhs = [1 K 1.7*K -(K+1)];            % given polynomial (1)
    assert(max(abs(lhs-rhs)) < 1e-12, 'Reconstruction mismatch at K=%.2f', K)
end
disp('Closed-loop polynomial reconstruction verified for sample K values.')

% Response at the main-band edge theta = pi  (z = e^{j*pi} = -1)
G_edge = polyval(num_GOL,-1) / polyval(den_GOL,-1);
fprintf('G_OL(-1) = %.4f  (%.2f dB)\n', G_edge, 20*log10(abs(G_edge)))

% Resonant digital frequency of the complex-conjugate pole pair
theta_res = angle(pp(abs(imag(pp))>1e-9));
fprintf('Resonant digital frequency = %.4f rad = %.1f deg\n', ...
    abs(theta_res(1)), abs(theta_res(1))*180/pi)
