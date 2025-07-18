function display_amode(ax, x, y, y_lim)
    delete(findobj(ax, 'Tag', 'plot_amode'));
    plot(ax, x, y, '-r', 'LineWidth', 2, 'Tag', 'plot_amode');
    axis(ax, 'tight');
    ylim(ax, [0 y_lim]);
end

