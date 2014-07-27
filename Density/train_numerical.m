N = 512; % sparse vector length
T = 5;  % sparsity
K = 512; % signal length

[x y_true y sigma A mask Num_obs sigma_true] = signal_generator_exp(N, T, K);
Num_tot = zeros(size(Num_obs));

[y_BP x_BP t_BP] = my_density_csf('BP', x, y, T, sigma, A, mask, 0.6); %epsilon(d,1));
%[y_DBP x_DBP t_DBP] = my_density_csf('DBP', x, y, T, sigma, A, mask, 15, 0, Num_obs, Num_tot, 0);
% [y_TBP x_TBP t_TBP] = my_density_csf('TBP', x, y, T, sigma, A, mask, epsilon(d,3), 0.065, Num_obs, Num_tot);
% [y_BCS x_BCS t_BCS] = my_density_csf('BCS', x, y, T, sigma, A, mask);
% [y_BOUND x_BOUND t_BOUND] = my_density_csf('DBP', x, y, T, sigma, A, mask, epsilon(d,4) , 0, Num_obs, Num_tot, C(d), sigma_true);
        
Ey_BP = norm(y_true-y_BP)/norm(y_true)
% Ey_BCS = norm(y_true-y_BCS)/norm(y_true)
% Ey_DBP = norm(y_true-y_DBP)/norm(y_true)      
% Ey_TBP = norm(y_true-y_TBP)/norm(y_true)
% Ey_BOUND = norm(y_true-y_BOUND)/norm(y_true)

figure
plot(x);
hold on
plot(x_BP, 'r');
