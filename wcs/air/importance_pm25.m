clear
clc

pm25i = load('.\airdata\pm25i.txt');
fprintf('load finished\n');

dim = size(pm25i);
n = dim(2);
snapn = dim(1);
is_train = 1;
train_per = 0.7;

ks = 10:10:50;
ms = 40:40:200;
if (is_train == 0) base = base_dct4(n);
else bases = train_pm25(pm25i', train_per, ks); 
end
fprintf('training finishedn\n');
    


% cost = load('d:\document\codes\data\cost\cost.stationu64.txt');
% cost = load('d:\document\codes\data\cost\cost.gradient.32.32.txt');
% cost = load('.\airdata\3gcost.32.32.txt');
% cost = cost(1:32, 1:32);
% costv = hilbertcurve(cost);
% alp = 1;
% prob = costv.^(-alp)/sum(costv.^(-alp));
prob = ones(1, 1024) / 1024;

ures = zeros(length(ms), 2);
bres = zeros(length(ms), 2);
sampling = 2;
denoise = 0;
alpha = 1;
times = 10;
snap = 1;

for j = 33:33  % which day
    for i = 1 : length(ms)  % number of samples
        m = ms(i);
        k = ks(i);
        usnr = 0;
        bsnr = 0;
        for t = 1 : times;  % number of iterations
            usamplex = randperm(n, m);
            bsamplex = [];
            prob = ones(1, 1024) ./ 1024;
            mm_step = 20;
            for mm = mm_step: mm_step: m
                kk = mm / 4;
                mm_add = mm_step;
                bsamplex_add = randsamplewtr(n, mm_add, prob);
                bsamplex = [bsamplex, bsamplex_add]; 
                if (is_train == 1) 
                    base = bases(i,:,:);
                end
                recover_b = my_csf(pm25i(j, bsamplex)', bsamplex, base, kk , 'OMP', denoise);
                recover_b(find(recover_b < 0)) = 0;
                total_imp = sum(recover_b .^ (alpha) ) - sum(recover_b(bsamplex) .^ (alpha));
                prob = recover_b' .^(alpha) ./ total_imp;
                prob(bsamplex) = 0;
            end
            if (is_train == 1)
                base = bases(length(ms), :, :);
            end
            recover_u = my_csf(pm25i(j, usamplex)', usamplex, base, k, 'OMP', denoise);
            recover_b = my_csf(pm25i(j, bsamplex)', bsamplex, base, k, 'OMP', denoise);
            usnr = usnr + exp_snr_imp(pm25i(j,:)', recover_u, 0);
            bsnr = bsnr + exp_snr_imp(pm25i(j,:)', recover_b, 0);
        end
        ures(i, :) = [m usnr/times/snap];
        bres(i, :) = [m bsnr/times/snap];
        fprintf('%d\n', m);
    end
end
figure
plot(ures(:,1),ures(:,2));
hold on
plot(bres(:,1),bres(:,2),'r');
