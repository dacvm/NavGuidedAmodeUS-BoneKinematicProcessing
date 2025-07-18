function display_amode_window(ax, lb, m, ub)
    delete(findobj(ax, 'Tag', 'plot_amodewindow'));
    xline(ax, m, '--b', 'LineWidth', 1, 'Tag', 'plot_amodewindow');
    xline(ax, lb, '-b', 'LineWidth', 1, 'Tag', 'plot_amodewindow');
    xline(ax, ub, '-b', 'LineWidth', 1, 'Tag', 'plot_amodewindow');
end

