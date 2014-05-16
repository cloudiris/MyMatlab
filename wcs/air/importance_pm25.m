% clear
% clc
%  
% pm25i = load('.\airdata\pm25i.txt');
% fprintf('load finished\n');
% 
% pos = randsample(1024, 300);
% pm25i_train = pm25i(1:500, 1:300);
% pm25i = pm25i(:,1:300);
% dim = size(pm25i);
% n = dim(2);
% snapn = dim(1);
% is_train = 1;
% km_ratio = 5;
% 
% ms = 40:20:160;
% ks = ms / km_ratio;
% stairs = zeros(length(ms), 1);
% m_stair = zeros(length(ms), 10);
% for i = 1:length(ms)
%     m_stair(i, 1) = 40;
%     j = 1;
%     tot = 40;
%     while (tot < ms(i))
%         j = j + 1;
%         m_stair(i, j) = 20;
%         tot = tot + 20;
%     end
%     stairs(i) = j;
% end
% 
% 
% 
% if (is_train == 0) base = base_dct4(n);
% else bases = train_pm25(pm25i_train', ks); 
% end
% fprintf('training finishedn\n');
    

prob = ones(1, n) / n;

ures = zeros(length(ms), 2);
ures_normal = zeros(length(ms), 2);
bres = zeros(length(ms), 2);
bres_normal = zeros(length(ms), 2);

sampling = 2;
denoise = 1;
alpha = 2;
times = 30;
snap = 1;

for j = 30:30 % which day
    for i = 1 : length(ms)  % number of samples
        m = ms(i);
        k = ks(i);
        usnr = 0;
        bsnr = 0;
        usnr_normal = 0;
        bsnr_normal = 0;
        for t = 1 : times;  % number of iterations
            bsamplex = [];
            prob = ones(1, n) ./ n;
            kk = 0;
            for sr = 1:stairs(i) 
                mm_add = m_stair(i, sr);
                kk = kk + mm_add / km_ratio;
                bsamplex_add = randsamplewtr(n, mm_add, prob);
                bsamplex = [bsamplex, bsamplex_add]; 
                if (is_train == 1) 
                     base = reshape(bases(:,:,i), size(bases(:,:,i),1), size(bases(:,:,i),2));
                end
                recover_b = my_csf(pm25i(j, bsamplex)', bsamplex, base, kk , 'OMP', denoise);
                recover_b(find(recover_b < 0)) = 0;
                importance = recover_b;
                total_imp = sum(importance .^ (alpha)) - sum(importance(bsamplex) .^ (alpha));
                prob = importance' .^(alpha) ./ total_imp;
                prob(bsamplex) = 0;
            end
            if (is_train == 1)
                 base = reshape(bases(:,:,length(ms)), size(bases(:,:,length(ms)),1), size(bases(:,:,length(ms)),2));
            end
            usamplex = randperm(n, m);
            recover_u = my_csf(pm25i(j, usamplex)', usamplex, base, k, 'OMP', denoise);
            recover_b = my_csf(pm25i(j, bsamplex)', bsamplex, base, k, 'OMP', denoise);
            usnr = usnr + exp_snr_imp(pm25i(j,:)', recover_u, 0);
            bsnr = bsnr + exp_snr_imp(pm25i(j,:)', recover_b, 0);
            %usnr_normal = usnr_normal + sum(abs(recover_u - pm25i(j, :)')) / n;
            %bsnr_normal = bsnr_normal + sum(abs(recover_b - pm25i(j, :)')) / n;
        end
        ures(i, :) = [m usnr/times/snap];
        bres(i, :) = [m bsnr/times/snap];
%         ures_normal(i, :) = [m usnr_normal/times/snap];
%         bres_normal(i, :) = [m bsnr_normal/times/snap];
        fprintf('%d\n', m);
    end
end
figure
plot(ures(:,1),ures(:,2));
hold on
plot(bres(:,1),bres(:,2),'r');
% hold on 
% plot(ures_normal(:,1),ures_normal(:,2),'g');
% hold on
% plot(bres_normal(:,1),bres_normal(:,2),'y');
