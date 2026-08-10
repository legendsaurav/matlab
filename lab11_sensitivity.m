%% EXPERIMENT 11: FULL-RANGE SENSITIVITY ANALYSIS (-100% to +100%)
% Varies each K element linearly from 0 to 2*nominal (i.e., -100% to +100% change)
% NO EARLY BREAK — records all points, stable or unstable

clear; close all; clc;

%% SYSTEM MATRICES
A = [1.1,  0.3, -1.5;
     0.1,  3.5,  2.2;
     0.4,  2.4, -1.1];
B = [0.1,  0.0;
     0.0,  1.1;
     1.0,  1.0];
C = [1.0,  2.0,  0.0;
     0.0,  0.0,  1.0];

K_nominal = [ 0.195155,  -0.484232,  -1.826888;
             -15.859694,  -5.224297,  10.454099];

elem_names = {'K11','K12','K13','K21','K22','K23'};
outdir = pwd;
num_steps = 200;  % 200 points = 1% resolution per step

summary = cell(length(elem_names), 10);

for idx = 1:6
    if idx <= 3
        i = 1; j = idx;
    else
        i = 2; j = idx - 3;
    end
    
    elem_name = elem_names{idx};
    nom = K_nominal(i,j);
    
    fprintf('Processing %s (nominal = %.6f)...\n', elem_name);
    
    % Build linear sweep: from -100% to +100% of nominal
    % K_value = nom * (1 + pct/100), where pct goes from -100 to +100
    pct_range = linspace(-100, 100, num_steps);
    
    all_rows = cell(num_steps, 10);
    
    for s = 1:num_steps
        pct = pct_range(s);
        kval = nom * (1 + pct/100);
        
        K = K_nominal;
        K(i,j) = kval;
        Acl = A - B*K;
        ev = eig(Acl);
        maxmag = max(abs(ev));
        isStable = all(abs(ev) < 1);
        
        all_rows{s,1} = s;
        all_rows{s,2} = elem_name;
        all_rows{s,3} = kval;
        all_rows{s,4} = kval - nom;
        all_rows{s,5} = pct;
        all_rows{s,6} = ev(1);
        all_rows{s,7} = ev(2);
        all_rows{s,8} = ev(3);
        all_rows{s,9} = maxmag;
        all_rows{s,10} = ternary(isStable, 'YES', 'NO');
    end
    
    T = cell2table(all_rows, 'VariableNames', ...
        {'Step_Number','K_Element_Name','K_Value','Change_from_Nominal','Change_Percent',...
         'Eigenvalue_1','Eigenvalue_2','Eigenvalue_3','Max_Magnitude','Stable'});
    
    % Write CSV
    fname = fullfile(outdir, sprintf('sensitivity_FULL_%s.csv', elem_name));
    writetable(T, fname);
    fprintf('  Wrote %s (%d rows)\n', fname, num_steps);
    
    % Find stability boundaries from the full data
    stable_mask = strcmp(T.Stable, 'YES');
    stable_pct = pct_range(stable_mask);
    
    if ~isempty(stable_pct)
        lower_bound = nom * (1 + min(stable_pct)/100);
        upper_bound = nom * (1 + max(stable_pct)/100);
        total_width = upper_bound - lower_bound;
        max_neg = abs(min(stable_pct));
        max_pos = abs(max(stable_pct));
    else
        lower_bound = NaN; upper_bound = NaN; total_width = 0;
        max_neg = 0; max_pos = 0;
    end
    
    % Most sensitive eigenvalue at boundary
    if ~isempty(stable_pct)
        bound_idx = find(stable_mask, 1, 'last'); % upper boundary
        if ~isempty(bound_idx) && bound_idx < num_steps
            ev_unstable = T.Eigenvalue_3(bound_idx+1); % approximate
        else
            ev_unstable = T.Eigenvalue_1(bound_idx);
        end
    else
        ev_unstable = NaN;
    end
    
    denom = abs(nom); if denom < eps, denom = eps; end
    critical = ternary(total_width/denom < 0.10, 'YES', 'NO');
    
    summary{idx,1} = elem_name;
    summary{idx,2} = nom;
    summary{idx,3} = lower_bound;
    summary{idx,4} = upper_bound;
    summary{idx,5} = total_width;
    summary{idx,6} = max_neg;
    summary{idx,7} = max_pos;
    summary{idx,8} = ev_unstable;
    summary{idx,9} = critical;
    summary{idx,10} = sum(stable_mask);
    
    % ========== PLOT ==========
    fig = figure('Visible','off');
    
    % Color by stability
    stable_x = pct_range(stable_mask);
    stable_y = T.Max_Magnitude(stable_mask);
    unstable_x = pct_range(~stable_mask);
    unstable_y = T.Max_Magnitude(~stable_mask);
    
    plot(stable_x, stable_y, 'g-o', 'LineWidth', 1.2, 'MarkerSize', 3); hold on;
    plot(unstable_x, unstable_y, 'r-x', 'LineWidth', 1.2, 'MarkerSize', 3);
    
    yline(0.5, 'b--', 'LineWidth', 1.5);
    yline(1.0, 'm--', 'LineWidth', 1.5);
    xline(0, 'k:', 'LineWidth', 0.8);
    
    xlabel('Percentage change from nominal (%)');
    ylabel('Max |eigenvalue|');
    title(sprintf('%s: Full Sweep -100%% to +100%% (Stable: %.1f%% to +%.1f%%)', ...
        elem_name, -max_neg, max_pos));
    legend('Stable','Unstable','Nominal max|λ|=0.5','Stability boundary |λ|=1',...
        'Nominal K','Location','best');
    grid on;
    
    % Auto-scale y-axis but keep minimum at 0
    ymax = max(T.Max_Magnitude) * 1.05;
    if ymax < 2, ymax = 2; end
    ylim([0, ymax]);
    
    saveas(fig, fullfile(outdir, sprintf('plot_FULL_%s.png', elem_name)));
    close(fig);
end

%% MASTER SUMMARY
sum_tbl = cell2table(summary, 'VariableNames', ...
    {'K_Element','Nominal_Value','Lower_Boundary','Upper_Boundary',...
     'Stable_Range_Width','Max_Neg_Change_pct','Max_Pos_Change_pct',...
     'Most_Sensitive_Eigenvalue','Critical','Stable_Points_Count'});
writetable(sum_tbl, fullfile(outdir, 'sensitivity_FULL_master_summary.csv'));
fprintf('\nMaster summary: sensitivity_FULL_master_summary.csv\n');

%% RANKED BY SENSITIVITY
widths = cell2mat(summary(:,5));
[~, order] = sort(widths);
ranked = summary(order, :);
rank_tbl = cell2table(ranked, 'VariableNames', ...
    {'K_Element','Nominal_Value','Lower_Boundary','Upper_Boundary',...
     'Stable_Range_Width','Max_Neg_Change_pct','Max_Pos_Change_pct',...
     'Most_Sensitive_Eigenvalue','Critical','Stable_Points_Count'});
writetable(rank_tbl, fullfile(outdir, 'sensitivity_FULL_ranked.csv'));
fprintf('Ranked table: sensitivity_FULL_ranked.csv\n');

disp('COMPLETE: Full-range sensitivity analysis done.');

%% Helper
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end