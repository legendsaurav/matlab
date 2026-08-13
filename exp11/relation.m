T = readtable("D:\matlab\artifacts\pole_K_data_short.csv");

% --------------------------------------------------
% Independent variables
% --------------------------------------------------

a0 = T.a0;
a1 = T.a1;
a2 = T.a2;

% Sort poles so pole ordering does not affect the test
P = sort([T.p1 T.p2 T.p3], 2);

p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);

% --------------------------------------------------
% K columns
% --------------------------------------------------

K = T{:, startsWith(T.Properties.VariableNames, "k")};
Knames = T.Properties.VariableNames( ...
    startsWith(T.Properties.VariableNames, "k"));

% --------------------------------------------------
% Linear models
% --------------------------------------------------

% K = b0 + b1*a0 + b2*a1 + b3*a2
Xa = [ones(height(T),1) a0 a1 a2];

% K = b0 + b1*p1 + b2*p2 + b3*p3
Xp = [ones(height(T),1) p1 p2 p3];

% --------------------------------------------------
% Calculate R^2 for every K
% --------------------------------------------------

R2_coeff = zeros(1,size(K,2));
R2_poles  = zeros(1,size(K,2));

for i = 1:size(K,2)

    y = K(:,i);

    % ---- coefficients ----
    b = Xa \ y;
    ypred = Xa*b;

    R2_coeff(i) = ...
        1 - sum((y-ypred).^2) / ...
            sum((y-mean(y)).^2);

    % ---- poles ----
    b = Xp \ y;
    ypred = Xp*b;

    R2_poles(i) = ...
        1 - sum((y-ypred).^2) / ...
            sum((y-mean(y)).^2);
end

% --------------------------------------------------
% Display results
% --------------------------------------------------

Result = table(Knames', R2_coeff', R2_poles', ...
    'VariableNames', {'K','R2_coefficients','R2_poles'});

disp(Result)\
re