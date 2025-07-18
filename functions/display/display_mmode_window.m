function display_mmode_window(ax, lb, m, ub)
    delete(findobj(ax, 'Tag', 'plot_mmodewindow'));
    yline(ax, m, '--w', 'LineWidth', 1, 'Tag', 'plot_mmodewindow');
    yline(ax, lb, '-w', 'LineWidth', 1, 'Tag', 'plot_mmodewindow');
    yline(ax, ub, '-w', 'LineWidth', 1, 'Tag', 'plot_mmodewindow');
end

