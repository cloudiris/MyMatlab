rates = [0.1:0.1:1.0];
n_round = 5;
n_rate = size(rates,2);
n_method = 6;
bestthr = [15 15 15 15 15 15 15 15 15 15];
bestthr2 = [20 20 20 20 20 20 20 20 20 20];

e_y = zeros(n_rate, n_method);

for i_rate = 1:n_rate
    for i_round = 1:n_round
        [x y_true y sigma A mask Num_obs Num_tot] = signal_readin_traffic(rates(i_rate));

        N = size(x, 1);
        T = 10;
        K = size(y, 1);
        %e_BP = finding_epsilon('BP', 0, 1, 0.01, 10, x, y_true, y, sigma, A, mask, Num_obs, Num_tot, 0)
        [y_BP x_BP t_BP] = my_density_csf('BP', x, y, T, sigma, A, mask, 0.2);
        
        %e_DBP = finding_epsilon('DBP', 0, 0.1, 0.001, 10, x, y_true, y, sigma, A, mask, Num_obs, Num_tot, 0)
        [y_DBP x_DBP t_DBP] = my_density_csf('DBP', x, y, T, sigma, A, mask, 0.02, 0, Num_obs, Num_tot);
        
        %e_TBP = finding_epsilon('TBP', 0, 0.1, 0.001, 10, x, y_true, y, sigma, A, mask, Num_obs, Num_tot, bestthr(i_rate))
        [y_TBP x_TBP t_TBP] = my_density_csf('TBP', x, y, T, sigma, A, mask, 0.015, bestthr(i_rate), Num_obs, Num_tot);
        
        %e_TDBP = finding_epsilon('TDBP', 0, 0.1, 0.001, 10, x, y_true, y, sigma, A, mask, Num_obs, Num_tot, bestthr2(i_rate))
        [y_TDBP x_TDBP t_TDBP] = my_density_csf('TDBP', x, y, T, sigma, A, mask, 0.02, bestthr2(i_rate), Num_obs, Num_tot, y_true);
        
        % [y_OMP x_OMP t_OMP] = my_density_csf('OMP', x, y, T, sigma, A, mask);
        % [y_DOMP x_DOMP t_DOMP] = my_density_csf('DOMP', x, y, T, sigma, A, mask);
        % [y_TDOMP x_TDOMP t_TDOMP] = my_density_csf('TDOMP', x, y, T, sigma, A, mask, 0, 0.5);
        [y_BCS x_BCS t_BCS] = my_density_csf('BCS', x, y, T, sigma, A, mask);

        % reconstruction error
        E_BP = norm(x-x_BP)/norm(x);
        E_BCS = norm(x-x_BCS)/norm(x);
        E_DBP = norm(x-x_DBP)/norm(x);
        E_TDBP = norm(x-x_TDBP)/norm(x);
        E_TBP = norm(x-x_TBP)/norm(x);

        % E_OMP = norm(x-x_OMP)/norm(x);
        % E_DOMP = norm(x-x_DOMP)/norm(x);
        % E_TDOMP = norm(x-x_TDOMP)/norm(x);

        e_y(i_rate, 1) = e_y(i_rate, 1) + norm(y_true-y_BP)/norm(y_true) / n_round;
        e_y(i_rate, 2) = e_y(i_rate, 2) + norm(y_true-y_TBP)/norm(y_true) / n_round;
        e_y(i_rate, 3) = e_y(i_rate, 3) + norm(y_true-y_DBP)/norm(y_true) / n_round;
        e_y(i_rate, 4) = e_y(i_rate, 4) + norm(y_true-y_TDBP)/norm(y_true) / n_round;
        e_y(i_rate, 6) = e_y(i_rate, 6) + norm(y_true-y_BCS)/norm(y_true) / n_round;
        
        Ey_BP = norm(y_true-y_BP)/norm(y_true)
        Ey_BCS = norm(y_true-y_BCS)/norm(y_true)
        Ey_DBP = norm(y_true-y_DBP)/norm(y_true)
        Ey_TDBP = norm(y_true-y_TDBP)/norm(y_true)
        Ey_TBP = norm(y_true-y_TBP)/norm(y_true)

        % Ey_OMP = norm(y_true-y_OMP)/norm(y_true);
        % Ey_DOMP = norm(y_true-y_DOMP)/norm(y_true);
        % Ey_TDOMP = norm(y_true-y_TDOMP)/norm(y_true);
    end
end

figure;
plot(rates, e_y(:,1));
hold on
plot(rates, e_y(:,2), 'r');
hold on
plot(rates, e_y(:,3), 'g');
hold on
plot(rates, e_y(:,4), 'm');
hold on
plot(rates, e_y(:,6), 'y');