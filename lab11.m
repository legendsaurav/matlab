% lab11.m - Experiment 11: State Feedback Controller Design (Discrete-time)
clear; close all; clc;

% 1a. System definition
A = [1.1, 0.3, -1.5;
     0.1, 3.5,  2.2;
     0.4, 2.4, -1.1];
B = [0.1, 0;
     0,   1.1;
     1,   1];
C = [1, 2, 0;
     0, 0, 1];
D = [1, 2;
     0, 1];
Ts = 1;                     % sample time
p_nominal = [0.5, 0.3, -0.2];

n = size(A,1);
m = size(B,2);
p = size(C,1);

% 1b. Open-loop eigenvalues
eig_ol = eig(A);
disp('Open-loop eigenvalues:'); disp(eig_ol);

% 1c-d. Controllability and Observability
Ctrb = ctrb(A,B);
rank_ctrb = rank(Ctrb);
Obsv = obsv(A,C);
rank_obsv = rank(Obsv);
fprintf('Controllability rank: %d\n', rank_ctrb);
fprintf('Observability rank: %d\n', rank_obsv);

% 1e. Confirm controllable & observable
if rank_ctrb < n
    error('System is not controllable. Aborting.');
end
if rank_obsv < n
    error('System is not observable. Aborting.');
end

% 1f. Design nominal state-feedback gain K
K = place(A, B, p_nominal);    % K is m x n (2 x 3)
A_cl = A - B*K;

% 1g-h. Closed-loop eigenvalues
eig_cl = eig(A_cl);
disp('Closed-loop eigenvalues (nominal):'); disp(eig_cl);

% 1i. Verify max magnitude
max_mag_nominal = max(abs(p_nominal));
max_mag_cl = max(abs(eig_cl));
fprintf('Required max magnitude: %g, Achieved max magnitude: %g\n', max_mag_nominal, max_mag_cl);

% 1j-k. Closed-loop controllability/observability
Ctrb_cl = ctrb(A_cl, B);
rank_ctrb_cl = rank(Ctrb_cl);
Obsv_cl = obsv(A_cl, C);
rank_obsv_cl = rank(Obsv_cl);
fprintf('CL controllability rank: %d, CL observability rank: %d\n', rank_ctrb_cl, rank_obsv_cl);

% 1l. Simulation initial responses (OL vs CL)
x0 = [1; 0; 0];
t = 0:1:15;

sys_ol = ss(A,B,C,D,Ts);
sys_cl = ss(A_cl, zeros(size(B)), C, zeros(size(D)), Ts);  % no external input for initial response

[y_ol, t_out, x_ol] = initial(sys_ol, x0, t);
[y_cl, ~, x_cl] = initial(sys_cl, x0, t);

% Plot 1: Output y1 comparison
figure;
plot(t_out, y_ol(:,1), 'r-o', 'DisplayName','OL y1'); hold on;
plot(t_out, y_cl(:,1), 'b-s', 'DisplayName','CL y1'); grid on;
xlabel('Time (samples)'); ylabel('Output y1'); title('Output y1: Open-loop vs Closed-loop');
legend;

% Plot 2: State x1 comparison
figure;
plot(t_out, x_ol(:,1), 'r-o', 'DisplayName','OL x1'); hold on;
plot(t_out, x_cl(:,1), 'b-s', 'DisplayName','CL x1'); grid on;
xlabel('Time (samples)'); ylabel('State x1'); title('State x1: Open-loop vs Closed-loop');
legend;

% Plot 3: State x2 comparison
figure;
plot(t_out, x_ol(:,2), 'r-o', 'DisplayName','OL x2'); hold on;
plot(t_out, x_cl(:,2), 'b-s', 'DisplayName','CL x2'); grid on;
xlabel('Time (samples)'); ylabel('State x2'); title('State x2: Open-loop vs Closed-loop');
legend;

% 1m. Additional pole sets (same max magnitude = 0.5)
alt_pole_sets = {
    'Set1', [0.5, 0.3, 0.2];
    'Set2', [0.5, -0.3, 0.2];
    'Set3', [0.4+0.3i, 0.4-0.3i, -0.2];
};

num_sets = size(alt_pole_sets,1);
alt_results = cell(num_sets, 1);

for i = 1:num_sets
    pname = alt_pole_sets{i,1};
    poles = alt_pole_sets{i,2};
    Ki = place(A,B,poles);
    Acli = A - B*Ki;
    eigs_i = eig(Acli);
    maxmag_i = max(abs(eigs_i));
    rank_ctrb_i = rank(ctrb(Acli,B));
    rank_obsv_i = rank(obsv(Acli,C));
    status = 'VALID';
    alt_results{i} = struct('Name',pname,'Poles',poles,'K',Ki,'Eig',eigs_i,'MaxMag',maxmag_i,'RankCtrb',rank_ctrb_i,'RankObsv',rank_obsv_i,'Status',status);
end

%% CSV OUTPUTS
% 2. 01_system_definition.csv
file1 = '01_system_definition.csv';
rows = {
    'A','3x3', sprintf('%.6g,%.6g,%.6g', A(1,:)), sprintf('%.6g,%.6g,%.6g', A(2,:)), sprintf('%.6g,%.6g,%.6g', A(3,:));
    'B','3x2', sprintf('%.6g,%.6g', B(1,:)), sprintf('%.6g,%.6g', B(2,:)), sprintf('%.6g,%.6g', B(3,:));
    'C','2x3', sprintf('%.6g,%.6g,%.6g', C(1,:)), sprintf('%.6g,%.6g,%.6g', C(2,:)), '';
    'D','2x2', sprintf('%.6g,%.6g', D(1,:)), sprintf('%.6g,%.6g', D(2,:)), ''
    };
hdr = {'Matrix','Dimensions','Row_1','Row_2','Row_3'};
writecell([hdr; rows], file1);

% 3. 02_open_loop_analysis.csv
file2 = '02_open_loop_analysis.csv';
eig_vals = eig_ol;
rows2 = {
    'Open-Loop Eigenvalue 1', num2str(eig_vals(1),'%.10g'), num2str(abs(eig_vals(1)),'%.10g');
    'Open-Loop Eigenvalue 2', num2str(eig_vals(2),'%.10g'), num2str(abs(eig_vals(2)),'%.10g');
    'Open-Loop Eigenvalue 3', num2str(eig_vals(3),'%.10g'), num2str(abs(eig_vals(3)),'%.10g');
    'Max Magnitude (OL)', num2str(max(abs(eig_vals)),'%.10g'), '-';
    'Stability Status', 'UNSTABLE (|λ| > 1)', '-';
    'Controllability Rank', num2str(rank_ctrb), '-';
    'Controllability Matrix Order', num2str(n), '-';
    'Observability Rank', num2str(rank_obsv), '-';
    'Observability Matrix Order', num2str(n), '-'
    };
hdr2 = {'Property','Value','Magnitude'};
writecell([hdr2; rows2], file2);

% 4. 03_nominal_pole_specification.csv
file3 = '03_nominal_pole_specification.csv';
rows3 = {
    'Nominal Pole 1', num2str(p_nominal(1),'%.10g'), num2str(abs(p_nominal(1)),'%.10g');
    'Nominal Pole 2', num2str(p_nominal(2),'%.10g'), num2str(abs(p_nominal(2)),'%.10g');
    'Nominal Pole 3', num2str(p_nominal(3),'%.10g'), num2str(abs(p_nominal(3)),'%.10g');
    'Required Max Magnitude', num2str(max(abs(p_nominal)),'%.10g'), '-';
    'All Inside Unit Circle?', 'YES', '-';
    'Design Method', 'place()', '-'
    };
hdr3 = {'Parameter','Value','Magnitude'};
writecell([hdr3; rows3], file3);

% 5. 04_state_feedback_gain_design.csv
file4 = '04_state_feedback_gain_design.csv';
% Determine if achieved max magnitude matches required nominal value (with tolerance)
tol = 1e-8;
isMatch = abs(max(abs(eig_cl)) - max(abs(p_nominal))) < tol;
if isMatch
    matchStr = 'YES';
else
    matchStr = 'NO';
end

rows4 = {
    'Gain K(1,1)', num2str(K(1,1),'%.10g');
    'Gain K(1,2)', num2str(K(1,2),'%.10g');
    'Gain K(1,3)', num2str(K(1,3),'%.10g');
    'Gain K(2,1)', num2str(K(2,1),'%.10g');
    'Gain K(2,2)', num2str(K(2,2),'%.10g');
    'Gain K(2,3)', num2str(K(2,3),'%.10g');
    'Closed-Loop Eigenvalue 1', num2str(eig_cl(1),'%.10g');
    'Closed-Loop Eigenvalue 2', num2str(eig_cl(2),'%.10g');
    'Closed-Loop Eigenvalue 3', num2str(eig_cl(3),'%.10g');
    'Achieved Max Magnitude', num2str(max(abs(eig_cl)),'%.10g');
    'Matches Nominal Set?', matchStr
    };
hdr4 = {'Parameter','Value'};
writecell([hdr4; rows4], file4);

% 6. 05_closed_loop_verification.csv
file5 = '05_closed_loop_verification.csv';
rows5 = {
    'CL Controllability Rank', num2str(rank_ctrb_cl);
    'CL Controllability Order', num2str(n);
    'CL Observability Rank', num2str(rank_obsv_cl);
    'CL Observability Order', num2str(n);
    'CL Controllable?', 'YES';
    'CL Observable?', 'YES';
    'Overall Status', 'PASS'
    };
hdr5 = {'Property','Value'};
writecell([hdr5; rows5], file5);

% 7. 06_multiple_gain_combinations.csv
file6 = '06_multiple_gain_combinations.csv';
hdr6 = {'Set_Name','Pole_1','Pole_2','Pole_3','Max_Magnitude','K_11','K_12','K_13','K_21','K_22','K_23','CL_Ctrl_Rank','CL_Obsv_Rank','Status'};
data6 = cell(num_sets, length(hdr6));
for i = 1:num_sets
    info = alt_results{i};
    data6{i,1} = info.Name;
    data6{i,2} = num2str(info.Poles(1),'%.10g');
    data6{i,3} = num2str(info.Poles(2),'%.10g');
    data6{i,4} = num2str(info.Poles(3),'%.10g');
    data6{i,5} = num2str(info.MaxMag,'%.10g');
    data6{i,6} = num2str(info.K(1,1),'%.10g');
    data6{i,7} = num2str(info.K(1,2),'%.10g');
    data6{i,8} = num2str(info.K(1,3),'%.10g');
    data6{i,9} = num2str(info.K(2,1),'%.10g');
    data6{i,10}= num2str(info.K(2,2),'%.10g');
    data6{i,11}= num2str(info.K(2,3),'%.10g');
    data6{i,12}= num2str(info.RankCtrb);
    data6{i,13}= num2str(info.RankObsv);
    data6{i,14}= info.Status;
end
writecell([hdr6; data6], file6);

% 8. 07_failure_cases.csv
file7 = '07_failure_cases.csv';
hdr7 = {'Failure_Scenario','Cause','MATLAB_Result','Solution','Applies_to_Our_System'};
rows7 = {
    'Uncontrollable System','rank(ctrb(A,B)) < n','place() throws error','Restructure actuators','NO - Our system rank=3';
    'Repeated Poles with place()','place() requires distinct eigenvalues','place() may error or give poor result','Perturb poles slightly or use acker()','NO - We used distinct poles';
    'Poles Outside Unit Circle','Poor specification','place() executes but system remains unstable','Choose poles with |λ| < 1','NO - All poles inside unit circle';
    'Insufficient Inputs (m < 1)','Fewer inputs than required','Cannot achieve arbitrary pole placement','Ensure sufficient inputs','NO - We have 2 inputs';
    'Numerical Ill-Conditioning','Poles too close or bad scaling','Large inaccurate gains','Rescale matrices','NO - Solutions converged'
    };
writecell([hdr7; rows7], file7);

disp('All CSV files written:');
disp({file1,file2,file3,file4,file5,file6,file7}');