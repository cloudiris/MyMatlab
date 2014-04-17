% clear
clear res* sres*
clc

%% pretreatment
% load data
% pm25i = load('pm25i.txt');
dim = size(pm25i);
n = dim(2);
d = sqrt(n);
baksnaps = load('sparsity_pm25.txt');
snapn = 1;
snapstart = randi(size(pm25i, 1));
% snapstart = 3666;
% snaps = sort(randsample(dim(1), snapn));
snaps = snapstart:snapstart + snapn - 1;
% snaps = baksnaps(snapstart:snapstart + snapn - 1);
% snaps = [1 3666];
fprintf('data loaded: %dx%d\n', dim(1), dim(2));
fprintf('snapshots: ');
disp(snaps);

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
m_max = n;
m_gap = 100;
ms = m_min : m_gap : m_max;
mn = length(ms);
b = base_dct4(n);
samplings = [1 2 3]; %0 ca, 1 random, 2 pw, 3 greedy, 4 llc
samplingn = length(samplings);
denoise = 0;
tolerance = 25;
times = 30; % single snapshot repeating
totalloopn = samplingn * mn * times * snapn;

alp = 0.1;

% to store

res = zeros(totalloopn, 5); % m|cost%|acc%|sampling|snap
res_random = zeros(mn, 3); % m|cost%|acc%
res_group = res_random;
res_ca = res_random;
res_greedy = res_random;
res_llc = res_random;

%% main loop
fprintf('main loop started\n');
loopcount = 0;

for samplingi = 1 : samplingn;
    sampling = samplings(samplingi);
    
    for mi = 1 : mn
        m = ms(mi);
        k = max(1, floor(m/log2(n)));
%         k = 50;
        
        for repeati = 1 : times
            %% sampling
            samplex = zeros(1, m);
            
            if (sampling == 0) %ca
                samplex = randsamplewtr(n, m, costv.^(-alp)/sum(costv.^(-alp)));
                
            elseif (sampling == 1) %random
                samplex = randperm(n, m);
                
            elseif (sampling == 2) %group
                m0 = m;
                if (2*m0 > n) m = n - m; end
                samplex = zeros(1, m);
                cands = randperm(n, 2*m);
                cox = mod(cands, 32);
                coy = ceil(cands./32);
                [tx1 tx2] = meshgrid(cox, cox);
                [ty1 ty2] = meshgrid(coy, coy);
                dis = (tx2-tx1).^2 + (ty2-ty1).^2 + diag(inf*ones(1, 2*m));
                for sam = 1 : m
                    [mins,mxs] = min(dis);
                    [~,myi] = min(mins);
                    mxi = mxs(myi);
                    idx1 = hcindex(32, cox(mxi), coy(mxi));
                    idx2 = hcindex(32, cox(myi), coy(myi));
                    if (costv(idx1) < costv(idx2))
                        if (2*m0 > n) samplex(sam) = idx2;
                        else samplex(sam) = idx1; end
                    else
                        if (2*m0 > n) samplex(sam) = idx1;
                        else samplex(sam) = idx2; end
                    end
                    dis([mxi myi], :) = inf*ones(2, 2*m);
                    dis(:, [mxi myi]) = inf*ones(2*m, 2);
                end
                if (2*m0 > n)
                    bif = ones(1,n);
                    bif(samplex) = 0;
                    samplex = find(bif==1);
                end
                m = m0;
                
            elseif (sampling == 3) %greedy
                [~, in] = sort(costv);
                samplex = in(1:m);
                
            elseif (sampling == 4) %llc
                bak_costv = costv;
                seeds = sort(randperm(n, m));
                l = 7;
                ll = ceil(l/2) - 1;
                rl = l - ll - 1;
                for sam = 1 : m
                    left = max(1, seeds(sam) - ll);
                    right = min(n, seeds(sam) + rl);
                    [~, offset] = min(bak_costv(left:right));
                    samplex(sam) = offset + left - 1;
                    bak_costv(samplex(sam)) = inf;
                end
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
                res(loopcount, :) = [m, costratio, acc, sampling, snap];
            end
%             acc = acc / snapn;
            
        end
        fprintf('sampling(%d)|m(%d)|k(%d) is over\n', sampling, m, k);
    end
    fprintf('sampling %d is over\n', sampling);
end

for mi = 1 : mn
    m = ms(mi);
    res_ca(mi, :) = mean(res(res(:,1)==m & res(:,4)==0, 1:3));
    res_random(mi, :) = mean(res(res(:,1)==m & res(:,4)==1, 1:3));
    res_group(mi, :) = mean(res(res(:,1)==m & res(:,4)==2, 1:3));
    res_greedy(mi, :) = mean(res(res(:,1)==m & res(:,4)==3, 1:3));
    res_llc(mi, :) = mean(res(res(:,1)==m & res(:,4)==4, 1:3));
end
if (length(ms) == size(res_greedy, 1) && length(ms) == size(res_random, 1))
    res_optimal = [ms' res_greedy(:,2) res_random(:,3)]; % random acc + greedy cost
end

figure
hold on
if (numel(res_random)>0) plot(res_random(:,2), res_random(:,3)); end
if (numel(res_group)>0) plot(res_group(:,2), res_group(:,3), 'r'); end
if (numel(res_ca)>0) plot(res_ca(:,2), res_ca(:,3), 'm'); end
if (numel(res_greedy)>0) plot(res_greedy(:,2), res_greedy(:,3), 'k'); end
if (numel(res_llc)>0) plot(res_llc(:,2), res_llc(:,3), 'g'); end
if (exist('res_optimal','var')) plot(res_optimal(:,2), res_optimal(:,3),'y'); end

imethod = 'linear';
accstops = 0:0.01:1;
sres_target = zeros(length(accstops), 1);

for samplingi = 1 : samplingn
    sampling = samplings(samplingi);
    if (sampling == 0) tres = res_ca;
    elseif (sampling == 1) tres = res_random;
    elseif (sampling == 2) tres = res_group;
    elseif (sampling == 3) tres = res_greedy;
    elseif (sampling == 4) tres = res_llc;
    end
    [sacc, idx] = sort(tres(:,3));
    [usacc, uidx] = unique(sacc);
    idx = idx(uidx);
    sres_target = interp1([0 tres(idx,3)' 1], [0 tres(idx,2)' 1], accstops, imethod)';
    if (sampling == 0) sres_ca = sres_target;
    elseif (sampling == 1) sres_random = sres_target;
    elseif (sampling == 2) sres_group = sres_target;
    elseif (sampling == 3) sres_greedy = sres_target;
    elseif (sampling == 4) sres_llc = sres_target;
    end
end

% interpolate optimal
if (exist('res_optimal','var'))
    tres = res_optimal;
    [sacc, idx] = sort(tres(:,3));
    [usacc, uidx] = unique(sacc);
    idx = idx(uidx);
    sres_target = interp1([0 tres(idx,3)' 1], [0 tres(idx,2)' 1], accstops, imethod)';
    sres_optimal = sres_target;
end

figure
hold on
if (exist('sres_random','var')) plot(accstops,sres_random); end
if (exist('sres_group','var')) plot(accstops,sres_group, 'r'); end
if (exist('sres_ca','var')) plot(accstops,sres_ca, 'm'); end
if (exist('sres_greedy','var')) plot(accstops,sres_greedy, 'k'); end
if (exist('sres_llc','var')) plot(accstops,sres_llc, 'g'); end
if (exist('sres_optimal','var')) plot(accstops, sres_optimal, 'y'); end

sres_all = [];
if (exist('sres_random','var')) sres_all = [sres_all sres_random]; end
if (exist('sres_group','var')) sres_all = [sres_all sres_group]; end
if (exist('sres_ca','var')) sres_all = [sres_all sres_ca]; end
if (exist('sres_greedy','var')) sres_all = [sres_all sres_greedy]; end
if (exist('sres_llc','var')) sres_all = [sres_all sres_llc]; end
if (exist('sres_optimal','var')) sres_all = [sres_all sres_optimal]; end
idx = [76 81 86 91 96];
figure;
bar(sres_all(idx, :));