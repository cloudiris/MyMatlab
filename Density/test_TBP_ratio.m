N = 512; % sparse vector length
T = 5;  % sparsity
K = 512; % signal length

Cs = 0:0.01:0.1;

n_C= size(Cs, 2);
error = zeros(n_C, 1);
for round = 1:1
[x y_true y sigma A mask Num_obs] = signal_generator_normal(N, T, K);
Num_tot = zeros(size(Num_obs));
    for i_C = 1:n_C
        [y_TBP x_TBP t_TBP] = my_density_csf('TBP_ratio', x, y, T, sigma, A, mask, 0.35, 0.7, Num_obs, Num_tot, Cs(i_C));
        error(i_C) = error(i_C) + norm(y_TBP - y) / norm(y)/ 10;
    end
end
figure;
plot(Cs, error);