
[x y_true y sigma A mask Num_obs Num_tot] = signal_readin_traffic(0.1, 0, 1);

N = size(x, 1);
T = 10;
K = size(y, 1);

Cs = 0:1000: 6000;
n_C = size(Cs, 2);
error = zeros(n_C, 1);

Max_CC = 10000;
CC = zeros(size(y(mask)));
s = sigma(mask);
d = Num_obs(mask);
for i = 1:size(y(mask), 1)
    if (d(i) > 1)     
        CC(i) = s(i) ^ 2 * (d(i) .^ 1.5) * ((d(i) .^ 2 - d(i)) / chi2inv(0.95, d(i) - 1) - 1);
    else
        CC(i) = Max_CC;
    end
end

figure;
hist(CC(i));

for i_C = 1:n_C
    [y_DBP x_DBP t_DBP] = my_density_csf('DBP', x, y, T, sigma, A, mask, 0.015, 0, Num_obs, Num_tot, Cs(i_C));
    error(i_C) = norm(y_DBP - y_true) / norm(y_true)
end

figure;
plot(Cs, error);




    