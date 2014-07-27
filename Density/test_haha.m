N = 512; % sparse vector length
T = 5;  % sparsity
K = 512; % signal length

[x y_true y sigma A mask Num_obs sigma_true] = signal_generator_normal(N, T, K);


[y_WBP x_WBP] = my_wbp_easy(y(mask), sigma(mask), A(mask, :), 5);

figure
plot(x);
hold on
plot(x_WBP, 'r');