function acc = recovery_throw_acc(y_observe, y, base, m, k, theta, guess, thr, denoise)
    acc = 0;
    n = size(y, 1);
    samplex = randperm(n, m);
    guess_vector = guess * ones(length(samplex), 1);
    r_vector = rand(length(samplex), 1);
    samplex_var = samplex(union(find(abs(y_observe(samplex) - y(samplex)) <= y(samplex) * thr), find(r_vector(:) > guess_vector(:)) ))';
%     samplex_var = [];
%     for i = 1:length(samplex)
%         if (abs(y(samplex(i)) - y_sample(i)) <= y(samplex(i)) * thr)
%             samplex_var = [samplex_var samplex(i)];
%         else
%             r = rand();
%             if (r >= guess) 
%                  samplex_var = [samplex_var samplex(i)];
%             end
%         end
%     end
    fprintf('simplex_var size\n');
    m
    k
    theta
    guess
    size(samplex_var)
    y_sample_var = y_observe(samplex_var);
    y_recover_var = my_csf(y_sample_var, samplex_var, base, k, 'OMP', denoise);
    %snr(i_m, i_theta) = snr(i_m, i_theta) + exp_snr(y, y_recover, 2); 
    %snr_var(i_m, i_theta) = snr_var(i_m, i_theta) + exp_snr(y, y_recover_var, 2);
    acc = exp_snr(y, y_recover_var, 2);