% ---------------------------------------------------------
% k11 - Six views in one figure
% ---------------------------------------------------------

figure('Color','white','Position',[100 100 1400 1000]);

% Color limits
limits = prctile(k11,[1 99]);

for i = 1:size(views,1)

    subplot(2,3,i);

    scatter3(a2, a1, a0, 8, k11, 'filled');

    xlabel('a_2');
    ylabel('a_1');
    zlabel('a_0');

    title(sprintf('(%c) View %d', 'a'+i-1, i));

    grid on;
    box on;

    clim(limits);
    colormap(turbo(256));

    view(views(i,1),views(i,2));

end

% One colorbar for the whole figure
cb = colorbar;
cb.Layout.Tile = 'east';

cb.Label.String = 'k_{11}';

sgtitle('Pole Space \rightarrow k_{11}', ...
        'FontSize',18, ...
        'FontWeight','bold');

% Save
exportgraphics(gcf,'k11_all_views.png','Resolution',300);

close(gcf);