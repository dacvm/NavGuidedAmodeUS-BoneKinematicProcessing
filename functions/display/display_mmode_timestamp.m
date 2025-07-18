function display_mmode_timestamp(ax, selected_timestamp)
    delete(findobj(ax, 'Tag', 'plot_mmodetimestamp'));
    xline(ax, selected_timestamp, 'r', 'LineWidth', 2, 'Tag', 'plot_mmodetimestamp');
end

