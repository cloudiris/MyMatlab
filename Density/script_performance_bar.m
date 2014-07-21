N = 512; % sparse vector length
T = 5;  % sparsity
K = 512; % signal length

error = zeros(3, 4);
n_round = 1;

epsilon = [0.66 7 0.7;
           0.66 7 0.7
           0.55 2 0.35];
for d = 1:3
    for i = 1:n_round

        if (d == 1) [x y_true y sigma A mask Num_obs] = signal_generator_normal(N, T, K);
        elseif (d == 2) [x y_true y sigma A mask Num_obs] = signal_generator_uniform(N, T, K);
        elseif (d == 3) [x y_true y sigma A mask Num_obs] = signal_generator_uniform(N, T, K);
        end    
        [y_BP x_BP t_BP] = my_density_csf('BP', x, y, T, sigma, A, mask, epsilon(d,1));
        [y_DBP x_DBP t_DBP] = my_density_csf('DBP', x, y, T, sigma, A, mask, epsilon(d,2), 0.5, Num_obs, Num_tot);
        [y_TBP x_TBP t_TBP] = my_density_csf('TBP', x, y, T, sigma, A, mask, epsilon(d,3), 2, Num_obs, Num_tot);
        [y_BCS x_BCS t_BCS] = my_density_csf('BCS', x, y, T, sigma, A, mask);

        Ey_BP = norm(y_true-y_BP)/norm(y_true)
        Ey_BCS = norm(y_true-y_BCS)/norm(y_true)
        Ey_DBP = norm(y_true-y_DBP)/norm(y_true)      
        Ey_TBP = norm(y_true-y_TBP)/norm(y_true)

        error(d,1) = error(d,1) + Ey_BP / n_round;
        error(d,2) = error(d,2) + Ey_BCS / n_round;
        error(d,3) = error(d,3) + Ey_DBP / n_round;
        error(d,4) = error(d,4) + Ey_TBP / n_round;
    end
end

%error(2,:) = error(1,:);
error(3,:) = error(1,:);

figure;
bar(error, 1);
legend('BP','BCS','WBP','TBP');
set(gca,'XTickLabel',{'Normal','Uniform','Exponential'});
set(gcf,'color','white')
applyhatch(gcf,'\.x.');
