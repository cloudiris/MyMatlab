%function [x y y_hat sigma_hat A mask] = signal_process(N, T, K)

function [x y_true y sigma A mask] = signal_process(y_true, Num_obs, Max_sigma, Max_y)
rand('state', 1);
randn('state', 2);
%
%N sparse vector length
%T sparsity
%K length of the original signal
%

mask = find(Num_obs(:) > 0);
y_hat = zeros(K, 1);
sigma_hat = zeros(K, 1);

for i = 1:K
    z = Z(i, :)';
    z = z(1: Num_obs(i));
    if (Num_obs(i) > 0) 
        y_hat(i) = mean(z);
        if (Num_obs(i) > 1)
            sigma_hat(i) = sqrt(Max_sigma / (Num_obs(i) ^ 2) + (mean(z .^ 2) - y_hat(i) .^ 2) / (Num_obs(i) - 1));
        else 
            sigma_hat(i) = Max_sigma;
        end
    else
        y_hat(i) = Max_y;
        sigma_hat(i) = Max_sigma;
    end 
end


