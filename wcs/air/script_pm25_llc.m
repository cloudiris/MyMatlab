% clear
% clc

%% pretreatment
% load data
% pm25i = load('pm25i.txt');
dim = size(pm25i);
n = dim(2);
d = sqrt(n);
baksnaps = load('sparsity_pm25.txt');
snapn = 1;
snapstart = 3666;
% snaps = sort(randsample(dim(1), snapn));
snaps = snapstart:snapstart + snapn - 1;
% snaps = baksnaps(snapstart:snapstart + snapn - 1);
% snaps = [1 3666];
fprintf('data loaded: %dx%d\n', dim(1), dim(2));

% load cost
% cost = load('d:\document\codes\data\cost\cost.stationu64.txt');
% cost = load('d:\document\codes\data\cost\cost.gradient.32.32.txt');
% cost = load('cost.gradient.exp.32.32.txt');
cost = load('D:\Document\Codes\Matlab\wcs\3gcost.32.32.txt');
% cost = load('cost.step.32.32.txt');
cost = cost(1:d, 1:d);
costv = hilbertcurve(cost);
sumcost = sum(costv);
fprintf('cost loaded\n');

% setting
% k = 50;
m_min = 50;
m_max = 50;
m_gap = 100;
ms = m_min : m_gap : m_max;
mn = length(ms);
b = base_dct4(n);
denoise = 0;
tolerance = 25;
times = 100; % single snapshot repeating
totaltimes = times*mn;

res = zeros(totaltimes, 4); % m|cost%|acc%|snap

%% main loop
fprintf('main loop started\n');
loopcount = 0;
    
for mi = 1 : mn
    m = ms(mi);
    k = max(1, floor(m/log2(n)));
%         k = 50;

    for repeati = 1 : times
        %% sampling
        samplex = zeros(1, m);
        
        bak_costv = costv;
        seeds = sort(randperm(n, m));
        l = 100;
        ll = ceil(l/2) - 1;
        rl = l - ll - 1;
        for sam = 1 : m
            left = max(1, seeds(sam) - ll);
            right = min(n, seeds(sam) + rl);
            [~, offset] = min(bak_costv(left:right));
            samplex(sam) = offset + left - 1;
            bak_costv(samplex(sam)) = inf;
        end
        
        costratio = sum(costv(samplex)) / sumcost;

        %% recovery
        for snapi = 1 : snapn
            loopcount = loopcount + 1;
            snap = snaps(snapi);
            s = pm25i(snap, :)';
            re = my_csf(s(samplex), samplex, b, k, 'OMP', denoise);
            re(re < 0) = 0;
            acc = sum(abs(re-s) <= tolerance) / n;
            res(loopcount, :) = [m, costratio, acc, snap];
        end
%             acc = acc / snapn;

    end
    fprintf('m(%d)|k(%d) is over\n', m, k);
end

for mi = 1 : mn
    m = ms(mi);
    res_llc(mi, :) = mean(res(res(:,1)==m, 1:3));
end

imethod = 'linear';
accstops = 0:0.01:1;
tres = res_llc;
[sacc, idx] = sort(tres(:,3));
[usacc, uidx] = unique(sacc);
idx = idx(uidx);
sres_llc = interp1([0 tres(idx,3)' 1], [0 tres(idx,2)' 1], accstops, imethod)';