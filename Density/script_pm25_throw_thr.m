% clear
% clc
%  
% pm25i = load('/Users/xiaohong/Documents/Research/Project_Git/Data/pm25i.txt');
fprintf('load finished\n');

dim = size(pm25i);
n = dim(2);
snapn = dim(1);

ms = 200:20:400;
mk_ratio = 4;

thetas = 0.3:0.1:0.3; % Noise level
guesses = 1.0:-0.1:0.0; %Predict accuracy 
guess = 1.0;
thrs = 0.1:0.2:0.5; %throw threshold 
denoise = 0;
snr = zeros(length(ms), length(thetas));
snr_var = zeros(length(ms), length(thetas), length(guesses));
base = base_dct4(n);
tot_iteration = 20;

for j = 30:30
    y = pm25i(j,:)';
    
    for i_m = 1:size(ms,2)
        m = ms(i_m);
        k = m / mk_ratio;
        
        for i_theta = 1:length(thetas)
            theta = rand(1,n) * thetas(i_theta);
            y_observe = gaussian_disturb_prob(y, theta);
            
            for i_thr = 1:length(thrs)
                thr = thrs(i_thr);
                for iteration = 1:tot_iteration
                    snr_var(i_m, i_theta, i_thr) = snr_var(i_m, i_theta, i_thr) + recovery_throw_acc_prob(y_observe, y, base, m, k, theta, guess, thr, denoise);
        
                    samplex = randperm(n, m);
                    y_sample = y_observe(samplex);
                    y_recover = my_csf(y_sample, samplex, base, k , 'OMP', denoise);
                    snr(i_m, i_theta) = snr(i_m, i_theta) + exp_snr(y, y_recover, 2); 
                    
                end
            end
        end
    end
end
snr = snr ./ tot_iteration ./ length(thrs);
snr_var = snr_var ./ tot_iteration;
figure
plot(ms,snr(:,1),'r');
hold on
plot(ms,snr_var(:,1,1),'g');
hold on
plot(ms, snr_var(:,1, 2), 'm');
hold on
plot(ms, snr_var(:,1, 3), 'b');
% hold on
% plot(ms, snr_var(:,1, 4), 'r*-');
% hold on
% plot(ms, snr_var(:,1, 5), 'g*-');
% hold on
% plot(ms, snr_var(:,1, 6), 'y*-');
% hold on
% plot(ms, snr_var(:,1, 7), 'b*-');
% hold on
% plot(ms, snr_var(:,1, 8), 'ro-');
% hold on
% plot(ms, snr_var(:,1, 9), 'go-');
% hold on
% plot(ms, snr_var(:,1, 10), 'bo-');
% hold on
% plot(ms, snr_var(:,1, 11), 'yo-');


