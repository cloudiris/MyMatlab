clear
clc
 
pm25i = load('.\airdata\pm25i.txt');
fprintf('load finished\n');

dim = size(pm25i);
n = dim(2);
snapn = dim(1);

ms = 100:20:400;
mk_ratio = 4;
thetas = 0.3:0.3;
denoise = 0;
thr = 0.2;
snr = zeros(length(ms), length(thetas));
snr_var = zeros(length(ms), length(thetas));
base = base_dct4(n);

for j = 30:30
    y = pm25i(j,:)';
    for i_theta = 1:length(thetas)
        theta = thetas(i_theta);
        for i_m = 1:size(ms,2)
            m = ms(i_m);
            k = m / mk_ratio;
            for iteration = 1:20 
                samplex = randperm(n, m);
                y_observe = gaussian_disturb(y, theta);
                y_sample = y_observe(samplex);
                error = abs(y(samplex) - y_sample);
                samplex_var = samplex(find(abs(y(samplex) - y_sample) < y(samplex) * thr))';
                y_sample_var = y_observe(samplex_var);
                y_recover = my_csf(y_sample, samplex, base, k , 'OMP', denoise);
                y_recover_var = my_csf(y_sample_var, samplex_var, base, k, 'OMP', denoise);
                snr(i_m, i_theta) = snr(i_m, i_theta) + exp_snr(y, y_recover, 2); 
                snr_var(i_m, i_theta) = snr_var(i_m, i_theta) + exp_snr(y, y_recover_var, 2);
            end
            snr(i_m, i_theta) = snr(i_m, i_theta) / 20;
            snr_var(i_m, i_theta) = snr_var(i_m, i_theta) / 20;
        end
    end
end

figure
plot(ms,snr(:,1),'r');
hold on
plot(ms,snr_var(:,1),'b');


