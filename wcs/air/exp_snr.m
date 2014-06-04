function snr = exp_snr(y, y_recover, mode)

mask = find(y(:) >= quantile(y(:), 0.0));
diff = abs(y - y_recover);

if (mode == 1)
    snr = sum(diff) / sum(y);
elseif (mode == 0)
    snr = numel(find(diff(mask) <= 10)) / size(diff(mask), 1);
elseif (mode == 2)
    snr = sum(diff .^ 2) / sum(y .^ 2);
end