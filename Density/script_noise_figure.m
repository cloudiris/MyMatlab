
N = 512; % sparse vector length
T = 5;  % sparsity
K = 512; % signal length

noises = 0:0.005:0.06;
n_noise = size(noises, 2);
error = zeros(n_noise, 5);
l2 = zeros(n_noise, 1);
n_round = 5;

Mean_obs = 100;
Num_obs = max(floor((randn(K, 1) + 3) * Mean_obs / 3), 0);

epsilon_BP = [0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7];
epsilon_DBP = [7 7 7 7 7 7 7 7 7 7 7 7 7];
C_DBP = [75 75 75 75 75 75 75 75 75 75 75 75 75];
epsilon_TBP = [0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35];
thr_TBP = [0.025 0.03 0.035 0.04 0.045 0.055 0.065 0.075 0.085 0.095 0.1 0.105 0.11];
epsilon_Bound = [12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5];


epsilon_BP = [0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7];
epsilon_DBP = [15 15 15 15 15 15 15 15 15 15 15 15 15];
C_DBP = [0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6];
epsilon_TBP = [0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35];
thr_TBP = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
epsilon_Bound = [12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5];


for d = 1:n_round
    
    co = zeros(T, 1);
    for i = 1:T
        co(i) = max(randn() + 1, 0);
    end
    x_base = zeros(N,1);
    q = randperm(N);
    x_base(q(1:T)) = sign(randn(T,1)); 
    
    A = randn(K,N);
    A = A./repmat(sqrt(sum(A.^2,2)),[1,N]);	
 
    
    for i_noise = 1:13
        noise = randn(N, 1) * noises(i_noise);
        [x y_true y sigma mask sigma_true] = signal_generator_normal_noise(x_base, A, N, T, K, noise, Num_obs);
        Num_tot = zeros(size(Num_obs));
        l2(i_noise) = norm(noise) / norm(x);

         [y_BP x_BP t_BP] = my_density_csf('BP', x, y, T, sigma, A, mask, epsilon_BP(i_noise)); 
         [y_DBP x_DBP t_DBP] = my_density_csf('DBP', x, y, T, sigma, A, mask, epsilon_DBP(i_noise), 0, Num_obs, Num_tot, C_DBP(i_noise));
         [y_TBP x_TBP t_TBP] = my_density_csf('TBP_ratio', x, y, T, sigma, A, mask, epsilon_TBP(i_noise), thr_TBP(i_noise), Num_obs, Num_tot);
         [y_BCS x_BCS t_BCS] = my_density_csf('BCS', x, y, T, sigma, A, mask);
%         [y_BOUND x_BOUND t_BOUND] = my_density_csf('DBP', x, y, T, sigma, A, mask, epsilon_Bound(i_noise), 0, Num_obs, Num_tot, 11, sigma_true);

%         figure
%         plot(x);
%         hold on;
%         plot(x_DBP, 'r');
%         
%         figure
%         plot(x);
%         hold on;
%         plot(x_BP, 'g');
        Ey_BP = norm(y_true-y_BP) / 512;
        Ey_BCS = norm(y_true-y_BCS) / 512;
        Ey_DBP = norm(y_true-y_DBP) / 512; 
        Ey_TBP = norm(y_true - y_TBP) / 512;
%        Ey_BOUND = norm(y_true-y_BOUND) / 512;

        error(i_noise,1) = error(i_noise,1) + Ey_BP / n_round;
        error(i_noise,2) = error(i_noise,2) + Ey_BCS / n_round;
        error(i_noise,3) = error(i_noise,3) + Ey_DBP / n_round;
        error(i_noise,4) = error(i_noise,4) + Ey_TBP / n_round;
%        error(i_noise,5) = error(i_noise,5) + Ey_BOUND / n_round;
    end
end

figure;
plot(l2, error(:,1));
hold on;
plot(l2, error(:,2), 'g');
hold on;
plot(l2, error(:,3), 'r');
hold on;
plot(l2, error(:,4), 'm');
% hold on;
% plot(l2, error(:,5), 'y');

