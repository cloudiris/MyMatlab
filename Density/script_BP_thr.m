
N = 512; % sparse vector length
T = 5;  % sparsity
K = 512; % signal length


epsilon = [0.66 7 0.7;
           0.66 7 0.7
           0.55 2 0.35];
       
thrs = 0:0.005:0.2;
n_thr = size(thrs,2)       
n_round = 1;
er_TBP = zeros(n_thr, 3);

for i = 1:n_round       
    [x y_true y sigma A mask Num_obs] = signal_generator_normal(N, T, K);
    for i_thr = 1:n_thr
        [y_TBP x_TBP t_TBP] = my_density_csf('TBP', x, y, T, sigma, A, mask, 0.35, thrs(i_thr), Num_obs, Num_tot);
        Ey_TBP = norm(y_true-y_TBP)/norm(y_true);
        er_TBP(i_thr, 1) = er_TBP(i_thr, 1) + Ey_TBP / n_round;
        i_thr
        Ey_TBP
    end
%     [x y_true y sigma A mask Num_obs] = signal_generator_uniform(N, T, K);
%     for i_thr = 1:n_thr
%         [y_TBP x_TBP t_TBP] = my_density_csf('TBP', x, y, T, sigma, A, mask, 0.35, thrs(i_thr), Num_obs, Num_tot);
%         Ey_TBP = norm(y_true-y_TBP)/norm(y_true);
%         er_TBP(i_thr, 2) = er_TBP(i_thr, 2) + Ey_TBP / n_round;
%         i_thr
%         Ey_TBP
%     end
%     [x y_true y sigma A mask Num_obs] = signal_generator_exp(N, T, K);
%     for i_thr = 1:n_thr
%         [y_TBP x_TBP t_TBP] = my_density_csf('TBP', x, y, T, sigma, A, mask, 0.35, thrs(i_thr), Num_obs, Num_tot);
%         Ey_TBP = norm(y_true-y_TBP)/norm(y_true);
%         er_TBP(i_thr, 3) = er_TBP(i_thr, 3) + Ey_TBP / n_round;
%         i_thr
%         Ey_TBP
%     end
%     
end

figure
plot(thrs, er_TBP(:,1));
hold on
plot(thrs, er_TBP(:,2),'m');
hold on
plot(thrs, er_TBP(:,3),'r');