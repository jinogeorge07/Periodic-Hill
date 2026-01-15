% --- Compute mask midline from the actual 1D mask_y being plotted ---
mask_norm = mask_y / max(mask_y);          % normalize to [0,1] for this x
sgn = mask_norm - 0.5;
idx = find(sgn(1:end-1) .* sgn(2:end) <= 0, 1, 'first');  % first crossing of 0.5

if ~isempty(idx)
    y1 = data_y(idx);      y2 = data_y(idx+1);
    m1 = mask_norm(idx);   m2 = mask_norm(idx+1);
    mask_y_midline = y1 + (0.5 - m1) * (y2 - y1) / (m2 - m1);  % linear interpolation
else
    mask_y_midline = NaN;  % no crossing found
end