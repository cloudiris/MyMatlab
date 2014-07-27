function [x y y_hat sigma_hat A mask Num_obs sigma] = signal_generator_uniform(N, T, K)

%rand('state', 1);
%randn('state', 2);
%
%N sparse vector length
%T sparsity
%K length of the original signal
%
% random +/- 1 signal

%Constant setting
Mean_obs = 100;
Max_sigma = 1.5;
No_y = 10;


co = zeros(T, 1);
for i = 1:T
    co(i) = max(randn() + 1, 0);
end
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = sign(randn(T,1)); 

%power-law sparsity
%x(q(1:T)) = co;

%not perfect sparse vector
%x = randn(N,1) .* 0.03;
%x(q(1:T)) = ones(T, 1) * 1;

% Dictionary = base_dct4(N);
% d_pos = sort(randsample(N, K));
% A = Dictionary(d_pos, :);

A = randn(K,N);
A = A./repmat(sqrt(sum(A.^2,2)),[1,N]);	
y = A*x;

%Sigma Setting
sigma = rand(K, 1) / 2;

%Observations
Num_obs = max(floor(randn(K, 1) * Mean_obs * 2), 0);

sigma_true = sigma ./ sqrt(Num_obs);

mask = find(Num_obs(:) > 0);
y_hat = zeros(K, 1);
sigma_hat = zeros(K, 1);

for i = 1:K
    z = randn(Num_obs(i), 1) * sigma(i,1) + y(i);
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
% 
% sigma = sigma_hat(mask);
% sigma_pen = sigma + ones(size(sigma)) * max(sigma) ./ (sqrt(Num_obs(mask)));
% figure;
% yy = abs(y_hat - y);
% yy = yy(mask);
% bar(yy);
% hold on
% plot(sigma_pen, 'r*');
% 
% figure
% hist(Num_obs);





