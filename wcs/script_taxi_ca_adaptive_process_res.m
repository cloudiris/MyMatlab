interpm = 'pchip';

acctointerp = 0.85 : 0.05 : 0.90;
alphatointerp = 0:0.3:100;
res_para = zeros(intvnn, numel(acctointerp));
res_para_rate = res_para;
res_cost = res_para;
res_alpha = zeros(numel(alphatointerp) + 1, numel(acctointerp));
res_alpha(size(res_alpha, 1), :) = acctointerp;
res_rate = zeros(numel(alphatointerp), numel(acctointerp));
for intvni = 1 : intvnn
    intvn = intvns(intvni);
    res_intvn = res(find(res(:,9)==intvn), 1:8);
    
    for alphai = 1 : alphan
        alpha = alphas(alphai);
        res_tmp = res_intvn(find(res_intvn(:,1) == alpha), 2:8);
        x = res_tmp(:,2);
        x = [x; 1];
        y = res_tmp(:,6);
        y = [y; 1];
        [x, in] = unique(x);
        y = y(in);
        icost = interp1(x, y, acctointerp, interpm);
        res_alpha(alphai, :) = icost;
        
        y = res_tmp(:,1);
        y = [y; 1];
        y = y(in);
        irate = interp1(x, y, acctointerp, interpm);
        res_rate(alphai, :) = irate;
    end    
    
    for i = 1 : numel(acctointerp)
        clear icost;
        icost = interp1(alphas, res_alpha(1:alphan, i), alphatointerp, interpm);
        res_alpha(1:numel(alphatointerp), i) = icost;
        [mv, in] = min(res_alpha(1:numel(alphatointerp), i));
        res_para(intvni, i) = alphatointerp(in);
        res_cost(intvni, i) = mv;
        
        irate = interp1(alphas, res_rate(1:alphan, i), alphatointerp, interpm);
        res_rate(:, i) = irate;
        res_para_rate(intvni, i) = irate(in);
    end
    
%     figure;
%     grid on;
%     set(gca,'linewidth',2);
%     plot(alphatointerp, res_alpha(1:numel(alphatointerp), :)');
end

intvtointerp = 1:130;
res_ipara = zeros(numel(intvtointerp), numel(acctointerp));
res_ipara_rate = res_ipara;
res_icost = res_ipara;
for i = 1 : numel(acctointerp)
    res_ipara(:,i) = interp1(intvns', res_para(:,i), intvtointerp, interpm);
    res_ipara_rate(:,i) = interp1(intvns', res_para_rate(:,i), intvtointerp, interpm);
    res_icost(:,i) = interp1(intvns', res_cost(:,i), intvtointerp, interpm);
end

% figure;
% hold on;
% grid on;
% set(gca,'linewidth',2);
% plot(intvns, res_para(:,4));

figure;
hold on;
grid on;
set(gca,'linewidth',2);
plot(intvtointerp, res_ipara(:,2));


res_para = [intvns' res_para];
res_ipara = [intvtointerp' res_ipara];
save('res_taxi_ca_adaptive_para.txt','res_para', '-ascii');
save('res_taxi_ca_adaptive_ipara.txt','res_ipara', '-ascii');