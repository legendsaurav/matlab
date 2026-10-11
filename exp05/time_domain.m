s = tf("s");
G = 10 * ((s + 2.5)^2) / ((s + 12) * (s^2 + 0.12));

K_vec = linspace(0.1, 30, 800);

ess = nan(size(K_vec));
ts = nan(size(K_vec));
mp = nan(size(K_vec));
meetsAll = false(size(K_vec));

for i = 1:numel(K_vec)
    K = K_vec(i);
    T = feedback(K * G, 1);

    if isstable(T)
        yFinal = dcgain(T);
        ess(i) = abs(1 - yFinal) * 100;

        info = stepinfo(T, 'SettlingTimeThreshold', 0.05);
        ts(i) = info.SettlingTime;
        mp(i) = info.Overshoot;

        meetsAll(i) = ess(i) <= 5 && ts(i) <= 10 && mp(i) <= 2;
    end
end

idxEss = find(ess <= 5, 1, 'first');
idxTs = find(ts <= 10, 1, 'first');
idxMp = find(mp <= 2, 1, 'first');
idxFinal = find(meetsAll, 1, 'first');

if isempty(idxFinal)
    error('No K in the search range satisfies all three specifications.');
end

K_ess = K_vec(idxEss);
K_ts = K_vec(idxTs);
K_mp = K_vec(idxMp);
K_final = K_vec(idxFinal);

T_final = feedback(K_final * G, 1);

plotSpec(K_vec, ess, 5,  'Steady-State Error (%)', 'e_{ss} vs K', idxEss);
plotSpec(K_vec, ts, 10,  'Settling Time (s)',       'T_s vs K',    idxTs);
plotSpec(K_vec, mp, 2,   'Peak Overshoot (%)',      'M_p vs K',    idxMp);

figure('Color', 'w');
step(T_final);
grid on;
set(findall(gcf, 'Type', 'line'), 'LineWidth', 3);
title(sprintf('Closed-Loop Step Response at K = %.2f', K_final));
xlabel('Time (s)');
ylabel('Output');

figure('Color', 'w');
vals = [K_ess, K_ts, K_mp, K_final];
cats = categorical({'e_{ss}', 'T_s', 'M_p', 'Chosen'});
bar(cats, vals, 'FaceColor', [0.2 0.6 0.8], 'LineWidth', 3);
hold on;
yline(K_final, 'r:', 'LineWidth', 3);
scatter(4, K_final, 80, 'filled', 'MarkerFaceColor', 'r');
grid on;
ylabel('Gain K');
title('First Acceptable Gains and Final Selection');

fprintf('First acceptable K for e_ss: %.2f\n', K_ess);
fprintf('First acceptable K for T_s : %.2f\n', K_ts);
fprintf('First acceptable K for M_p  : %.2f\n', K_mp);
fprintf('Final selected K            : %.2f\n', K_final);

function plotSpec(K_vec, metric, req, yLabelText, plotTitleText, firstIdx)
    figure('Color', 'w');
    plot(K_vec, metric, 'LineWidth', 3);
    hold on;
    yline(req, 'r:', 'LineWidth', 3);

    if ~isempty(firstIdx)
        kMin = K_vec(firstIdx);
        yMin = metric(firstIdx);

        plot(kMin, yMin, 'kx', 'MarkerSize', 16, 'LineWidth', 4);
        plot(kMin, yMin, 'o', 'MarkerSize', 10, 'LineWidth', 3, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y');

        yLimits = ylim;
        yOffset = 0.04 * range(yLimits);
        labelY = max(req, yMin) + yOffset;

        text(kMin, labelY, sprintf('K = %.2f, y = %.2f', kMin, yMin), ...
            'FontWeight', 'bold', 'Color', 'k', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'BackgroundColor', 'w', ...
            'Margin', 2, ...
            'Clipping', 'off');
    end

    grid on;
    xlabel('K');
    ylabel(yLabelText);
    title(plotTitleText);
    set(findall(gcf, 'Type', 'line'), 'LineWidth', 3);
end
