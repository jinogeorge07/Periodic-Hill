function fill_walls_hillp(ax, x_val, y1_vals, y2_vals, y)
% Robust opaque fill of solid regions for two walls.
% It automatically decides which curve is top and bottom at each x.

    % Decide vertical extent of the domain
    if nargin < 5 || isempty(y)
        % Use walls themselves to get min/max if y not provided
        y_bot = min([y1_vals(:); y2_vals(:)]);
        y_top = max([y1_vals(:); y2_vals(:)]);
    else
        y_bot = min(y(:));
        y_top = max(y(:));
    end

    % Ensure column vectors
    x_val   = x_val(:);
    y1_vals = y1_vals(:);
    y2_vals = y2_vals(:);

    % Per-x top/bottom (handles any “flipped” inputs)
    y_top_curve = max(y1_vals, y2_vals);
    y_bot_curve = min(y1_vals, y2_vals);

    % Slight endpoint padding to avoid microscopic gaps
    xpadL = x_val(1);
    xpadR = x_val(end);

    % ---- TOP SOLID: area above the upper wall ----
    patch(ax, ...
        [xpadL; x_val; xpadR; flipud(x_val)], ...
        [y_top;  y_top_curve; y_top;  y_top*ones(size(x_val))], ...
        1, 'FaceColor',[0.2 0.2 0.2], 'EdgeColor','none', 'Clipping','on');

    % ---- BOTTOM SOLID: area below the lower wall ----
    patch(ax, ...
        [xpadL; x_val; xpadR; flipud(x_val)], ...
        [y_bot;  y_bot_curve; y_bot;  y_bot*ones(size(x_val))], ...
        1, 'FaceColor',[0.2 0.2 0.2], 'EdgeColor','none', 'Clipping','on');
end
