function snr = exp_snr_imp(y, y_recover, mode)

importance = (y .^ 2) ./ sum(y .^ 2);
mask = find(y(:) >= quantile(y(:), 0.7));
diff = abs(y - y_recover);

if (mode == 1)
    snr = (importance)' * (diff);
elseif (mode == 0)
    snr = numel(find(diff(mask) <= 10)) / size(diff(mask), 1);
elseif (mode == 2)
    snr = (importance)' * (diff .^ 2);
end

