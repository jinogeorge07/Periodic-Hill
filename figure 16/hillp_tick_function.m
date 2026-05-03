if hasBWR, colormap(bluewhitered); end

    % geometry markers
     %plot(x_solid,y_solid,'.','Color','k');
%     if haveWalls
%         plot(x_val, y1_vals, 'k');
%         plot(x_val, y2_vals, 'k');
%     end

    xlabel('x','FontSize',32);
    ylabel('y','FontSize',32);

    % ----- colorbar with big labels -----
    c = colorbar;
    c.Ticks      = linspace(c.Limits(1), c.Limits(2), 6);
    c.TickLabels = compose('%.2f', c.Ticks);
    c.FontSize   = BIG_TICKS;

    % ============ FRAME + TICKS ============

    ax = gca;

    % basic axis styling
    ax.FontSize         = BIG_TICKS;
    ax.LineWidth        = 1.6;
    ax.TickLength       = [0.02 0.02];
    ax.XAxis.TickDirection = 'out';
    ax.YAxis.TickDirection = 'out';
    ax.XAxis.Exponent   = 0;
    ax.YAxis.Exponent   = 0;

    % X ticks: from xmin to Lx, NXT values
    xmin    = min(x(:));
    ax.XLim = [xmin Lx];
    ax.XTick = linspace(xmin, Lx, NXT);
    ax.XTickLabel = arrayfun(@(v) sprintf('%d', v), ax.XTick, ...
                             'UniformOutput', false);

    % Y ticks: NYT values between current limits
    ymin = min(y(:));
    ymax = max(y(:));
    ax.YLim  = [ymin ymax];
    ax.YTick = linspace(ymin, ymax, NYT);

    % shift Y tick *labels* by +1.5 (hillp offset)
    yt = ax.YTick;
    ax.YTickLabel = arrayfun(@(v) sprintf('%d', v + 1.5), yt, ...
                             'UniformOutput', false);

    % shrink and position axes + colorbar so labels fit in frame
    set(ax, 'Units','normalized', 'Position', AX_POS);
    set(c,  'Units','normalized', 'Position', CB_POS);

    % ======================================
