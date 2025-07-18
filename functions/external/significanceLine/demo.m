fig = figure('Units', 'centimeters');
fig.Position(3) = 3.8; % the width of the figure
fig.Position(4) = 3.8; % the width of the figure

ax = axes('Units', 'centimeters');
ax.Position(3) = 3; % the width of the axes
ax.Position(4) = 3; % the height of the axes

% set the data
x = [1, 2, 3, 4];
y = [0.6, 0.4, 0.6, 0.8];

% bar plot
bar(ax, x, y);

% significance line between 1 and 2
significanceLine(ax, 1, 2, 0.8,...
    'edgeLength', 0.05,...
    'fontSize', 18);

% significance line between 1 and 4
significanceLine(ax, 1, 4, 1,...
    'marker', '**',...
    'edgeLength', 0.05,...
    'fontSize', 18);

% set the limitations of the axis
xlim(ax, [0.5, 4.5]);
ylim(ax, [0, 1.1]);

