clear;
clc;
close all;

% =========================================================
% Component name
% =========================================================
comp_name = 'k21';

% =========================================================
% Image filenames
% =========================================================
files = cell(6,1);

for i = 1:6
    files{i} = sprintf('%s_view_%d.png', comp_name, i);
end

labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

% =========================================================
% Read images
% =========================================================
imgs = cell(6,1);

for i = 1:6
    imgs{i} = imread(files{i});
end

% =========================================================
% Make all images the same size
% =========================================================
H = min(cellfun(@(x) size(x,1), imgs));
W = min(cellfun(@(x) size(x,2), imgs));

for i = 1:6
    imgs{i} = imresize(imgs{i}, [H W]);
end

% =========================================================
% Create figure
% =========================================================
fig = figure( ...
    'Color','white', ...
    'Units','pixels', ...
    'Position',[100 100 1800 1200]);

t = tiledlayout(2,3);

% More space between panels
t.TileSpacing = 'loose';
t.Padding = 'compact';

% =========================================================
% Add images
% =========================================================
for i = 1:6

    ax = nexttile;

    imshow(imgs{i}, 'Parent', ax);

    axis(ax,'image');
    axis(ax,'off');

    % -----------------------------------------------------
    % Small label BELOW the image
    % -----------------------------------------------------
    text(ax, 0.5, -0.035, labels{i}, ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',8, ...
        'FontWeight','bold', ...
        'Clipping','off');

end

% =========================================================
% Save
% =========================================================
output_name = sprintf('%s_combined.png', comp_name);

exportgraphics(fig, ...
    output_name, ...
    'Resolution',300);

close(fig);

fprintf('Saved: %s\n', output_name);