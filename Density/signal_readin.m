function [x y_true y_hat sigma A mask] = signal_readin(sampling_rate)

y_true = load('/Users/xiaohong/Documents/Research/Project_Git/Data/s.txt');
Num_tot = load('/Users/xiaohong/Documents/Research/Project_Git/Data/c.txt');
obs = load('/Users/xiaohong/Documents/Research/Project_Git/Data/obs.txt');

y_hat = zeros(K, 1);
sigma_hat = zeros(K, 1);
Max_sigma = max(y_true);
Max_y = -1;

N_tot = sum(Num_tot);
N_sample = floor(N_tot * sampling_rate);
pos_sample = sort(randperm(N_tot, N_sample));
obs = obs(pos_sample,:);
Num_obs = zeros(K, 1);

for i = 1:K
    z = obs(find(obs(:,1) == i), 2);
    Num_obs(i) = size(z, 1);
    if (Num_obs(i) > 0) 
        y_hat(i) = mean(z);
        if (Num_obs(i) > 1)
            sigma_hat(i) = sqrt(Max_sigma / (Num_obs(i) ^ 4) + (mean(z .^ 2) - y_hat(i) .^ 2) / (Num_obs(i) - 1));
        else 
            sigma_hat(i) = Max_sigma;
        end
    else
        y_hat(i) = Max_y;
        sigma_hat(i) = Max_sigma;
    end 
end
mask = find(Num_obs(:) > 0);

figure;
plot(y_true);
hold on
plot(y_hat, 'r');
