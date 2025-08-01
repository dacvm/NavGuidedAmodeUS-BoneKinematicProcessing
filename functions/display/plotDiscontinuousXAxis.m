function plotDiscontinuousXAxis(ax, x1, y1, x2, y2, varargin)
% Two segments with a narrow visual gap, uniform ticks, and configurable gap shading.
% Name-Value:
%   GapFraction  - visible gap width as fraction of total span [0.05]
%   NumTicks     - total number of x-ticks across the axis [9]
%   LineWidth    - curve width [1.5]
%   Color        - RGB color (applied to both segments) []
%   LineStyle    - '-', '--', ':', '-.' ['-']
%   TickFormat   - sprintf format for labels (default '%.0f' => integer)
%   ShadeAlpha   - opacity of the gap shading, 0..1 [0.04]
%   ShadeColor   - RGB color of the gap shading [[0 0 0]]
%   ShowShade    - logical flag to draw the shading [true]

    p = inputParser;
    addParameter(p,'GapFraction',0.05,@(v)isnumeric(v)&&isscalar(v)&&v>=0);
    addParameter(p,'NumTicks',9,@(v)isnumeric(v)&&isscalar(v)&&v>=2);
    addParameter(p,'LineWidth',1.5,@(v)isnumeric(v)&&isscalar(v)&&v>0);
    addParameter(p,'Color',[],@(v)isnumeric(v)&&isvector(v)&&numel(v)==3);
    addParameter(p,'LineStyle','-',@(s)ischar(s)||isstring(s));
    addParameter(p,'TickFormat','%.0f',@(s)ischar(s)||isstring(s));
    addParameter(p,'ShadeAlpha',0.04,@(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
    addParameter(p,'ShadeColor',[0 0 0],@(v)isnumeric(v)&&isvector(v)&&numel(v)==3);
    addParameter(p,'ShowShade',true,@(v)islogical(v)||ismember(v,[0 1]));
    parse(p,varargin{:});
    gapFrac = p.Results.GapFraction; nTicks = p.Results.NumTicks;
    lw = p.Results.LineWidth; col = p.Results.Color;
    lstyle = char(p.Results.LineStyle); tfmt = char(p.Results.TickFormat);
    shadeA = p.Results.ShadeAlpha; shadeC = p.Results.ShadeColor;
    doShade = logical(p.Results.ShowShade) && shadeA > 0;

    % Ensure column vectors & order by x-range
    x1 = x1(:); y1 = y1(:); x2 = x2(:); y2 = y2(:);
    [a1,b1] = bounds(x1); [a2,b2] = bounds(x2);
    if a2 < a1
        [x1,x2] = deal(x2,x1); [y1,y2] = deal(y2,y1);
        [a1,b1,a2,b2] = deal(a2,b2,a1,b1);
    end

    hold(ax,'on');

    if a2 <= b1
        p1 = plot(ax,x1,y1,'LineWidth',lw,'LineStyle',lstyle);
        p2 = plot(ax,x2,y2,'LineWidth',lw,'LineStyle',lstyle);
        if ~isempty(col), set([p1 p2],'Color',col); end
        set(ax,'XTickLabelRotation',0); hold(ax,'off'); return
    end

    % Compress the gap
    L1 = b1 - a1; L2 = b2 - a2; realGap = a2 - b1;
    totSpan = L1 + realGap + L2;
    gap_vis = max(gapFrac*totSpan, eps);
    shift   = realGap - gap_vis;
    x2p     = x2 - shift;

    % Plot segments
    p1 = plot(ax,x1 ,y1,'LineWidth',lw,'LineStyle',lstyle);
    p2 = plot(ax,x2p,y2,'LineWidth',lw,'LineStyle',lstyle);
    if isempty(col), p2.Color = p1.Color; else, set([p1 p2],'Color',col); end

    % Limits in displayed coordinates
    xMin = a1; xL = b1; xR = a2 - shift; xMax = b2 - shift;
    xlim(ax,[xMin xMax]);

    % Uniform ticks across whole axis (skip inside the gap) with rounded labels
    ticksDisp = linspace(xMin, xMax, nTicks);
    ticksDisp = ticksDisp(~(ticksDisp > xL & ticksDisp < xR));
    labels = arrayfun(@(t) sprintf(tfmt, t + (t >= xR)*shift), ticksDisp, ...
                      'UniformOutput', false);
    set(ax,'XTick',ticksDisp,'XTickLabel',labels,'XTickLabelRotation',0);

    % Gap markers
    xline(ax,xL,'--','Color',[.35 .35 .35],'LineWidth',1);
    xline(ax,xR,'--','Color',[.35 .35 .35],'LineWidth',1);

    % --- Shading: extend far beyond, but clip to axes and restore ylim ---
    if doShade
        yl0   = ylim(ax);                 % original limits for the data
        span0 = max(1, diff(yl0));
        B     = 1e6 * span0;              % huge vertical margin
    
        h = patch(ax,[xL xR xR xL], [yl0(1)-B yl0(1)-B yl0(2)+B yl0(2)+B], ...
                  shadeC, ...
                  'FaceAlpha',shadeA, 'EdgeColor','none', ...
                  'HitTest','off','PickableParts','none', ...
                  'Clipping','on');       % <— keep it inside the axes box
    
        uistack(h,'bottom');              % keep grid/lines above the shade
        ylim(ax, yl0);                    % restore the original y-limits
    end

    hold(ax,'off');
end
