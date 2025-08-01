function fig = plotColoredTrace(A, B, ax, alphaVal, lineWidth, lineStyle)
% PLOTCOLOREDTRACE  Plot B vs. time, colouring line segments and markers
% by angle A, with adjustable transparency, line width & style.
%
%   fig = plotColoredTrace(A, B)
%   fig = plotColoredTrace(A, B, ax)
%   fig = plotColoredTrace(A, B, ax, alphaVal)
%   fig = plotColoredTrace(A, B, ax, alphaVal, lineWidth)
%   fig = plotColoredTrace(A, B, ax, alphaVal, lineWidth, lineStyle)
%
% Inputs (all optional after A,B):
%   ax         – axes handle (if omitted, a new figure+axes is created)
%   alphaVal   – scalar in [0,1] for both line & marker transparency
%   lineWidth  – thickness of each segment (default 1.5)
%   lineStyle  – e.g. '-', '--', ':', '-.'  (default '-')

    narginchk(2,6);
    
    %— create or validate axes
    if nargin<3 || isempty(ax) || ~isgraphics(ax,'axes')
        fig = figure;
        ax  = axes(fig);
    else
        fig = ancestor(ax,'figure');
    end
    hold(ax,'on');
    
    %— defaults
    if nargin<4 || isempty(alphaVal),  alphaVal  = 1;    end
    if nargin<5 || isempty(lineWidth), lineWidth = 1.5;  end
    if nargin<6 || isempty(lineStyle), lineStyle = '-';  end
    
    %— sanity checks
    assert(isvector(A) && isvector(B) && numel(A)==numel(B), ...
        'A and B must be same-length vectors.');
    assert(isnumeric(alphaVal) && alphaVal>=0 && alphaVal<=1, ...
        'alphaVal must be in [0,1].');
    
    %— build colormap index from A
    cmap = turbo(256);
    idx  = round( rescale(A,1,size(cmap,1)) );
    
    %— draw each coloured segment
    t = 1:numel(A);
    for k = 1:numel(A)-1
        h = plot(ax, t(k:k+1), B(k:k+1), ...
            'Color',      cmap(idx(k),:), ...
            'LineWidth',  lineWidth, ...
            'LineStyle',  lineStyle);
        h.Color(4) = alphaVal;  % set RGBA alpha
    end
    
    % %— overlay the markers, same alpha
    % sc = scatter(ax, t, B, 30, A, 'filled');
    % sc.MarkerFaceAlpha = alphaVal;
    % sc.MarkerEdgeAlpha = alphaVal;
    
    %— finish styling
    colormap(ax, cmap);
    caxis(ax, [min(A) max(A)]);
    cb = colorbar(ax);
    cb.Label.String = '← Extention | Angle (°) | Flexion →';

    grid(ax,'on');
    box(ax,'on');
    hold(ax,'off');
end
