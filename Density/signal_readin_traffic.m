function [x y_true y_hat sigma_hat A mask Num_obs Num_tot] = signal_readin_traffic(sampling_rate)

y_true = load('/Users/xiaohong/Documents/Research/Project_Git/Data/s.txt');
Num_tot = load('/Users/xiaohong/Documents/Research/Project_Git/Data/c.txt');
obs = load('/Users/xiaohong/Documents/Research/Project_Git/Data/obs.txt');

K = size(y_true, 1);
y_hat = zeros(K, 1);
sigma_hat = zeros(K, 1);
Max_sigma = (max(y_true) - min(y_true)) / 2;
No_y = -1;

% N_tot = sum(Num_tot);
% N_sample = floor(N_tot * sampling_rate);
% pos_sample = sort(randperm(N_tot, N_sample));
% obs = obs(pos_sample,:);
Num_obs = zeros(K, 1);

for i = 1:K
    z = obs(find(obs(:,1) == i), 2);
    Num_tot(i) = size(z, 1);
    Num_obs(i) = floor(rand() * Num_tot(i) * sampling_rate);
    if (Num_obs(i) < 0) Num_obs(i) = 0;
    elseif Num_obs(i) > Num_tot(i) Num_obs(i) = Num_tot(i);
    end
    z = z(randperm(Num_tot(i), Num_obs(i)));
    
    if (Num_obs(i) > 0) 
        y_hat(i) = mean(z);
        if (Num_obs(i) > 1)
            sigma_hat(i) = sqrt((mean(z .^ 2) - y_hat(i) .^ 2) / (Num_obs(i) - 1));
        else 
            sigma_hat(i) = Max_sigma;
        end
    else
        y_hat(i) = No_y;
        sigma_hat(i) = Max_sigma;
    end 
end
mask = find(Num_obs(:) > 0);

A = dctmtx(K);
x = A * y_true;
A = inv(A);

% sigma = sigma_hat(mask);
% sigma_pen = sigma + ones(size(sigma)) * max(sigma) ./ (sqrt(Num_obs(mask)));
% figure;
% yy = abs(y_hat - y_true);
% yy = yy(mask);
% bar(yy);
% hold on
% plot(sigma_pen, 'r*');
% 
% figure;
% hist(Num_obs);


