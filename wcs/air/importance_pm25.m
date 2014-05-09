clear
clc

pm25i = load('.\airdata\pm25i.txt');
dim = size(pm25i);
n = dim(2);
snapn = dim(1);

k = 50;
m = 200;
base = base_dct4(n);

% cost = load('d:\document\codes\data\cost\cost.stationu64.txt');
% cost = load('d:\document\codes\data\cost\cost.gradient.32.32.txt');
% cost = load('.\airdata\3gcost.32.32.txt');
% cost = cost(1:32, 1:32);
% costv = hilbertcurve(cost);
% alp = 1;
% prob = costv.^(-alp)/sum(costv.^(-alp));
prob = ones(1, 1024) / 1024;

ks = 10:10:100;
ms = 40:40:400;

ures = zeros(length(ms), 2);
bres = zeros(length(ms), 2);
sampling = 2;
denoise = 0;
alpha = 1;

fprintf('load finished\n');
for j = 30:30  % which day
    importance = pm25i(j, :) ./ sum(pm25i(j, :));
    for i = 1 : length(ms)  % number of samples
        m = ms(i);
        k = ks(i);
        usnr = 0;
        bsnr = 0;
        times = 20;
        snap = 1;
        for t = 1 : times;  % number of iterations
            usamplex = randperm(n, m);
            bsamplex = [];
            prob = ones(1, 1024) ./ 1024;
            for mm = 20:20:m
                kk = mm / 4;
                mm_add = 20;
                bsamplex_add = randsamplewtr(n, mm_add, prob);
                bsamplex = [bsamplex, bsamplex_add]; 
                recover_b = my_csf(pm25i(j, bsamplex)', bsamplex, base, kk , 'OMP', denoise);
                recover_b(find(recover_b < 0)) = 0;
                total_imp = sum(recover_b .^ (alpha) ) - sum(recover_b(bsamplex) .^ (alpha));
                prob = recover_b' .^(alpha) ./ total_imp;
                prob(bsamplex) = 0;
            end
            usnr = usnr + sum (importance' .* (abs(...
                my_csf(pm25i(j, usamplex)', usamplex, base, k, 'OMP', denoise)-...
                pm25i(j,:)'...
                ))/n );
            bsnr = bsnr + sum(importance' .* (abs(...
                my_csf(pm25i(j, bsamplex)', bsamplex, base, k, 'OMP', denoise)-...
                pm25i(j,:)'...
                ))/n);
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

% 
% 
% for i = 1 : length(ms)
%     m = ms(i);
%     usnr = 0;
%     bsnr = 0;
%     ucost = 0;
%     bcost = 0;
%     times = 10;
%     snap = 1;
%     
%     for j = 3 : 3
%         for t = 1 : times
%             % sampling
%             usamplex = randperm(n, m);
%             bsamplex = zeros(m ,1);
%             if (sampling == 1)
%             %% CA
%                 bsamplex = randsamplewtr(n, m, prob);
%             %% GROUP
%             elseif (sampling == 2)
%                 cands = randperm(n, 2*m);
%                 cox = mod(cands, 32);
%                 coy = ceil(cands./32);
%                 [tx1 tx2] = meshgrid(cox, cox);
%                 [ty1 ty2] = meshgrid(coy, coy);
%                 dis = (tx2-tx1).^2 + (ty2-ty1).^2 + diag(inf*ones(1, 2*m));
%                 for sam = 1 : m
%                     [mins,mxs] = min(dis);
%                     [~,myi] = min(mins);
%                     mxi = mxs(myi);
%                     idx1 = hcindex(32, cox(mxi), coy(mxi));
%                     idx2 = hcindex(32, cox(myi), coy(myi));
%                     if (costv(idx1) < costv(idx2))
%                         bsamplex(sam) = idx1;
%                     else
%                         bsamplex(sam) = idx2;
%                     end
%                     dis([mxi myi], :) = inf*ones(2, 2*m);
%                     dis(:, [mxi myi]) = inf*ones(2*m, 2);
%                 end
%             end
%             
%             usnr = usnr + numel(find(abs(...
%                 my_csf(pm25i(j, usamplex)', usamplex, b, k, 'OMP', denoise)-...
%                 pm25i(j,:)'...
%                 ) <= 50))/n;
%             bsnr = bsnr + numel(find(abs(...
%                 my_csf(pm25i(j, bsamplex)', bsamplex, b, k, 'OMP', denoise)-...
%                 pm25i(j,:)'...
%                 ) <= 50))/n;
%             ucost = ucost + sum(costv(usamplex))/sum(costv);
%             bcost = bcost + sum(costv(bsamplex))/sum(costv);
% %       usnr = usnr + snr(pm25i(j, :), my_csf(pm25i(j, usamplex)', usamplex, b, k, 'OMP', 0));
% %       bsnr = bsnr + snr(pm25i(j, :), my_csf(pm25i(j, bsamplex)', bsamplex, b, k, 'OMP', 0));
% %       fprintf('%d\n', j);
%         end
%     end
%     ures(i, :) = [ucost/times/snap usnr/times/snap];
%     bres(i, :) = [bcost/times/snap bsnr/times/snap];
%     fprintf('%d\n', m);
% end
% 
% figure
% plot(ures(:,1),ures(:,2));
% hold on
% plot(bres(:,1),bres(:,2),'r');