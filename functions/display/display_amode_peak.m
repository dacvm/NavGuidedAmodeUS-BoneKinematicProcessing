function amode_peak_display(ax, peak_mm, peak_val)
    % plot the peak
    delete(findobj(ax, 'Tag', 'plot_amodepeak'));
    plot(ax, peak_mm, peak_val, 'or', 'MarkerFaceColor', 'r', 'Tag', 'plot_amodepeak');
end

