%% beat CA!

clear
clc

s = load('beatca.s.txt');
n = length(s);
k = 5;
b = base_dct4(n);
cost = load('cost.gradient.exp.32.32.txt');
costv = cost(1:n);
sumcost = sum(costv);

ms = 10:10:500;
samplings = 0:2;
alp = 1;

times = 10;

totaltimes = length(ms) * length(samplings) * times;
res = zeros(totaltimes, 4); % m|cost|acc|sampling

loopcounter = 0;
for mi = 1 : length(ms)
    m = ms(mi);
    
    for samplingi = 1 : length(samplings)
        sampling = samplings(samplingi);
        
        for timei = 1 : times
            loopcounter = loopcounter + 1;
            
            if (sampling == 0)
                samplex = randperm(n, m);

            elseif (sampling == 1)
                samplex = randsamplewtr(n, m, costv.^(-alp)/sum(costv.^(-alp)));

            elseif (sampling == 2)
                group = sort(randsample(n, 2*m));
                for i = 1 : m
                    if (costv(group(2*i-1)) < costv(group(2*i)))
                        samplex(i) = group(2*i-1);
                    else
                        samplex(i) = group(2*i);
                    end
                end
            end
            
            costratio = sum(costv(samplex)) / sumcost;
            re = my_csf(s(samplex), samplex, b, k, 'OMP', 0);
            res(loopcounter, :) = [m, costratio, snr(s, re), sampling];
        end
    end
    fprintf('%d\n', m);
end

res_random = zeros(length(ms), 4);
res_ca = res_random;
res_group = res_random;
for mi = 1 : length(ms)
    m = ms(mi);
    res_random(mi, :) = mean(res(res(:,4) == 0 & res(:,1)==m, :));
    res_ca(mi, :) = mean(res(res(:,4) == 1 & res(:,1) == m, :));
    res_group(mi, :) = mean(res(res(:,4) == 2 & res(:,1) == m, :));
end
        
figure
hold on
plot(res_random(:,2),res_random(:,3));
plot(res_ca(:,2),res_ca(:,3), 'r');
plot(res_group(:,2),res_group(:,3), 'm');