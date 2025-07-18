function display_shadedError(timepoints, means, stds, varargin)
% DISPLAY_SHADEDERROR Creates a plot with mean lines and shaded standard deviation areas
%
% Usage:
%   plot_mean_with_shaded_std(timepoints, means, stds)
%   plot_mean_with_shaded_std(timepoints, means, stds, 'PropertyName', PropertyValue, ...)
%   plot_mean_with_shaded_std(ax, timepoints, means, stds, ...)
%
% Inputs:
%   ax         - (optional) axes object to plot on
%   timepoints - vector of time points (x-axis)
%   means      - matrix where each column is the mean for one condition
%   stds       - matrix where each column is the std for one condition
%
% Optional Name-Value Parameters:
%   'labels'     - cell array of condition labels (default: {'Condition 1', 'Condition 2', ...})
%   'colors'     - matrix of RGB colors [n_conditions x 3] or cell array of color specs
%   'alpha'      - transparency for shaded areas (default: 0.3)
%   'linewidth'  - width of mean lines (default: 2)
%   'xlabel'     - x-axis label (default: 'Time', ignored if axes object provided)
%   'ylabel'     - y-axis label (default: 'Signal', ignored if axes object provided)
%   'title'      - plot title (default: '', ignored if axes object provided)
%   'legend_loc' - legend location (default: 'northeast')
%
% Example:
%   % Generate sample data
%   t = 0:0.1:20;
%   mean1 = 0.2 * sin(t) + 0.1 * exp(-t/5);
%   std1 = 0.05 + 0.02 * abs(sin(2*t));
%   mean2 = 0.1 * cos(t/2) - 0.05 * t/20;
%   std2 = 0.03 + 0.01 * abs(cos(3*t));
%   
%   % Plot on new figure
%   plot_mean_with_shaded_std(t, [mean1' mean2'], [std1' std2'], ...
%       'labels', {'stim', 'cue'}, ...
%       'xlabel', 'timepoint', ...
%       'ylabel', 'signal');
%   
%   % Plot on existing axes
%   fig = figure;
%   ax = axes(fig);
%   xlabel(ax, 'My Custom X Label');
%   ylabel(ax, 'My Custom Y Label');
%   plot_mean_with_shaded_std(ax, t, [mean1' mean2'], [std1' std2'], ...
%       'labels', {'stim', 'cue'});

% Check if first argument is an axes object
if nargin >= 1 && isa(timepoints, 'matlab.graphics.axis.Axes')
    ax = timepoints;
    timepoints = means;
    means = stds;
    stds = varargin{1};
    varargin = varargin(2:end);
    use_existing_axes = true;
else
    ax = [];
    use_existing_axes = false;
end

% Parse input arguments
p = inputParser;
addRequired(p, 'timepoints', @isnumeric);
addRequired(p, 'means', @isnumeric);
addRequired(p, 'stds', @isnumeric);
addParameter(p, 'labels', {}, @iscell);
addParameter(p, 'colors', [], @(x) isnumeric(x) || iscell(x));
addParameter(p, 'alpha', 0.3, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'linewidth', 2, @(x) isnumeric(x) && x > 0);
addParameter(p, 'xlabel', 'Time', @ischar);
addParameter(p, 'ylabel', 'Signal', @ischar);
addParameter(p, 'title', '', @ischar);
addParameter(p, 'legend_loc', 'northeast', @ischar);

parse(p, timepoints, means, stds, varargin{:});

% Extract parsed values
timepoints = p.Results.timepoints(:);
means = p.Results.means;
stds = p.Results.stds;
labels = p.Results.labels;
colors = p.Results.colors;
alpha_val = p.Results.alpha;
linewidth = p.Results.linewidth;
xlabel_str = p.Results.xlabel;
ylabel_str = p.Results.ylabel;
title_str = p.Results.title;
legend_loc = p.Results.legend_loc;

% Validate dimensions
if size(means, 1) ~= length(timepoints)
    if size(means, 2) == length(timepoints)
        means = means';
        stds = stds';
    else
        error('Size mismatch between timepoints and means');
    end
end

if any(size(means) ~= size(stds))
    error('Means and stds must have the same dimensions');
end

n_conditions = size(means, 2);

% Set default labels
if isempty(labels)
    labels = cell(n_conditions, 1);
    for i = 1:n_conditions
        labels{i} = sprintf('Condition %d', i);
    end
end

% Set default colors (similar to the plot you showed)
if isempty(colors)
    if n_conditions <= 2
        colors = [0.2, 0.4, 0.8; 1.0, 0.6, 0.2]; % Blue and orange
    else
        colors = lines(n_conditions);
    end
elseif iscell(colors)
    % Convert cell array of color specs to numeric
    color_matrix = zeros(length(colors), 3);
    for i = 1:length(colors)
        if ischar(colors{i}) && length(colors{i}) == 1
            color_matrix(i, :) = rem(floor((strfind('kbgcrmyw', colors{i}) - 1) * [0.25 0.5 1]), 2);
        else
            color_matrix(i, :) = colors{i};
        end
    end
    colors = color_matrix;
end

% Create the plot
if ~use_existing_axes
    figure;
    ax = gca;
end

% Set axes as current
axes(ax);
hold(ax, 'on');

% Plot shaded areas and mean lines for each condition
handles = zeros(n_conditions, 1);
for i = 1:n_conditions
    mean_line = means(:, i);
    std_line = stds(:, i);
    
    % Calculate upper and lower bounds
    upper_bound = mean_line + std_line;
    lower_bound = mean_line - std_line;
    
    % Create shaded area
    x_fill = [timepoints; flipud(timepoints)];
    y_fill = [upper_bound; flipud(lower_bound)];
    
    % Remove NaN values for fill
    valid_idx = ~(isnan(x_fill) | isnan(y_fill));
    x_fill = x_fill(valid_idx);
    y_fill = y_fill(valid_idx);
    
    % Plot shaded area
    fill(ax, x_fill, y_fill, colors(i, :), 'FaceAlpha', alpha_val, ...
         'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % Plot mean line
    handles(i) = plot(ax, timepoints, mean_line, 'Color', colors(i, :), ...
                     'LineWidth', linewidth);
end

% Customize plot (only if not using existing axes)
if ~use_existing_axes
    xlabel(ax, xlabel_str);
    ylabel(ax, ylabel_str);
    if ~isempty(title_str)
        title(ax, title_str);
    end
end

% Add legend
legend(ax, handles, labels, 'Location', legend_loc);

% % Improve plot appearance
% grid(ax, 'on');
% grid(ax, 'minor');
% box(ax, 'on');

% Set axis properties for better appearance
ax.FontSize = 10;
ax.LineWidth = 1;

hold(ax, 'off');

end