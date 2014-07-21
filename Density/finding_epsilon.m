function epsilon = finding_epsilon(method, low, high, stop, T, x, y_true, y, sigma, A, mask, Num_obs, Num_tot, thr)

epsilon = high;
while (abs(high - low) > stop)
    mid = (high + low) / 2
    [yy xx tt] = my_density_csf(method, x, y, T, sigma, A, mask, mid, thr, Num_obs, Num_tot);
    xx_20 = xx(find(xx >= quantile(xx, 0.8)));
    energy_20 = xx_20' * xx_20
    energy = xx' * xx
    haha = norm(yy) / norm(y)
    if (energy_20 >= 0.8 * energy && norm(yy) / norm(y) < 2 && norm(yy) / norm(y) > 0.5) 
        epsilon = min(mid, epsilon);
        low = mid;
    else
        high = mid;
    end
end


