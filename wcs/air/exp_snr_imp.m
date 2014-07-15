function snr = exp_snr_imp(y, y_recover, mode)

%err = (max(y) - min(y)) / 50;
%((y - y_min) .^ 2) ./ sum((y - y_min) .^ 2);
importance = compute_impotrance(y);

mask = find(y(:) >= quantile(y(:), 0.5));
diff = abs(y - y_recover);

if (mode == 1)
    snr = sum(importance .* diff);
elseif (mode == 0)
    snr = numel(find(diff(mask) <= 2)) / size(diff(mask), 1);
elseif (mode == 2)
    snr = sum(importance .* diff .^ 2);
elseif (mode == 11)
    snr = sum(diff);
elseif (mode == 10)
    snr = numel(find(diff <= 5)) / size(diff, 1);
elseif (mode == 12)
    snr = sum(diff .^ 2);
end

