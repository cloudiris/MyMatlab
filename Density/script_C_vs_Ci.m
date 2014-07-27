
N = 512; % sparse vector length
T = 5;  % sparsity
K = 512; % signal length


% epsilon = [0.66 7 0.7;
%            0.66 7 0.7
%            0.55 2 0.35];
       
Cs = 0:100:1000;
n_C = size(Cs,2);      
n_round = 1;
er_DBP = zeros(n_C, 3);

for i = 1:n_round       
    [x y_true y sigma A mask Num_obs sigma_true] = signal_generator_normal(N, T, K);
    Num_tot = zeros(size(Num_obs));
    for i_C = 1:n_C
        [y_DBP x_DBP t_DBP] = my_density_csf('DBP', x, y, T, sigma, A, mask, 1, 0, Num_obs, Num_tot, Cs(i_C));
        Ey_DBP = norm(y_true-y_DBP)/norm(y_true);
        er_DBP(i_C, 1) = er_DBP(i_C, 1) + Ey_DBP / n_round;
        Cs(i_C)
        Ey_DBP
    end
%     [x y_true y sigma A mask Num_obs] = signal_generator_uniform(N, T, K);
%     for i_C = 1:n_C
%         [y_DBP x_DBP t_DBP] = my_density_csf('DBP', x, y, T, sigma, A, mask, 0.9, 0, Num_obs, Num_tot, Cs(i_C));
%         Ey_DBP = norm(y_true-y_DBP)/norm(y_true);
%         er_DBP(i_C, 2) = er_DBP(i_C, 2) + Ey_DBP / n_round;
%         Cs(i_C)
%         Ey_DBP
%     end
%     [x y_true y sigma A mask Num_obs] = signal_generator_exp(N, T, K);
%     for i_C = 1:n_C
%         [y_DBP x_DBP t_DBP] = my_density_csf('DBP', x, y, T, sigma, A, mask, 1, 0, Num_obs, Num_tot, Cs(i_C));
%         Ey_DBP = norm(y_true-y_DBP)/norm(y_true);
%         er_DBP(i_C, 3) = er_DBP(i_C, 3) + Ey_DBP / n_round;
%         Cs(i_C)
%         Ey_DBP
%     end
    
end

figure
plot(Cs, er_DBP(:,1));
hold on
plot(Cs, er_DBP(:,2),'m');
hold on
plot(Cs, er_DBP(:,3),'r');

CC = zeros(size(y(mask)));
s = sigma(mask);
d = Num_obs(mask);
for i = 1:size(y(mask), 1)
    CC(i) = s(i) ^ 2 * (d(i) .^ 1.5) * ((d(i) - 1) / chi2inv(0.05, d(i) - 1) - 1 / d(i));
end

figure;
hist(CC);
     