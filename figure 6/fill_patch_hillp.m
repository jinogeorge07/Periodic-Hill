% --- Fill region BELOW the hill y1_vals(x) with grey ---
ax = gca;
hold(ax,'on');

% bottom of the domain (0 in your plots)
y_bottom = ax.YLim(1);     % should be 0.0

% make sure column vectors
xcol = data_x(:);          % same x used for y1_vals
ycol = y1_vals(:);

% polygon: up along the hill, back along the bottom
x_fill = [xcol;      flipud(xcol)      ];
y_fill = [ycol;      y_bottom*ones(size(xcol))];

patch(x_fill, y_fill, [0.3 0.3 0.3], ...
      'EdgeColor','none', ...
      'FaceAlpha',1, ...
      'Parent',ax);