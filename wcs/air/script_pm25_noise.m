clear res res_* sres_*
clc

%% pretreatment
% load data
% pm25i = load('pm25i.txt');
dim = size(pm25i);
n = dim(2);
d = sqrt(n);
baksnaps = load('sparsity_pm25.txt');
snapn = 1;
snapstart = 3666;
snaps = sort(randsample(dim(1), snapn));
% snaps = snapstart:snapstart + snapn - 1;
% snaps = baksnaps(snapstart:snapstart + snapn - 1);
% snaps = [1 3666];
fprintf('data loaded: %dx%d\n', dim(1), dim(2));

% load cost
% cost = load('d:\document\codes\data\cost\cost.stationu64.txt');
% cost = load('d:\document\codes\data\cost\cost.gradient.32.32.txt');
% cost = load('cost.gradient.exp.32.32.txt');
cost = load('D:\Document\Codes\Matlab\wcs\3gcost.32.32.txt');
tmp = load('d:\document\codes\data\temp.grid.32.32.txt');
tmp = hilbertcurve(tmp);
% cost = load('cost.step.32.32.txt');
cost = cost(1:d, 1:d);
costv = hilbertcurve(cost);
sumcost = sum(costv);
fprintf('cost loaded\n');

noi = 10;

ntimes = 1;

alp = 0.6;
% prob = fcostv.^(-alp)/sum(fcostv.^(-alp));

% setting
k = 50;
m_min = 50;
m_max = n;
m_gap = 100;
ms = m_min : m_gap : m_max;
mn = length(ms);
b = base_dct4(n);
samplings = [1 2]; % 0 ca, 1 random, 2 pw, 3 greedy
samplingn = length(samplings);
denoise = 0;
tolerance = 25;
times = 10; % single snapshot repeating
totalloopn = samplingn * mn * times * snapn;

% to store

res = zeros(totalloopn, 6); % m|cost%|acc%|noise|sampling|snap
res_random = zeros(mn, 4); % m|cost%|acc%|noise
res_group = res_random;
res_ca = res_random;
res_greedy = res_random;

%% main loop
fprintf('main loop started\n');
loopcount = 0;

for ntimei = 1 : ntimes
    nois = noise_n(n, noi*norm(costv));
    fcostv = (nois'+costv);
    fcostv = rand(n, 1);
    prob = fcostv.^(-alp)/sum(fcostv.^(-alp));
    noisel = norm(nois)/norm(costv);

    for samplingi = 1 : samplingn;
        sampling = samplings(samplingi);

        for mi = 1 : mn
            m = ms(mi);
    %         k = max(1, floor(m/log2(n)));
    %         k = 50;

            for repeati = 1 : times
                %% sampling
    %             samplex = zeros(1, m);

                if (sampling == 0) %ca
                    samplex = randsamplewtr(n, m, prob);

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
                        if (randi(2)==1) samplex(sam) = idx1;
                        else samplex(sam) = idx2; end
%                         if (fcostv(idx1) < fcostv(idx2))
%                             if (2*m0 > n) samplex(sam) = idx2;
%                             else samplex(sam) = idx1; end
%                         else
%                             if (2*m0 > n) samplex(sam) = idx1;
%                             else samplex(sam) = idx2; end
%                         end
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
                    [~, in] = sort(fcostv);
                    samplex = in(1:m);

                end
                costratio = sum(costv(samplex)) / sumcost;

                %% recovery
                for snapi = 1 : snapn
                    loopcount = loopcount + 1;
                    snap = snaps(snapi);
                    s = pm25i(snap, :)';
                    s = tmp';
                    re = my_csf(s(samplex), samplex, b, k, 'OMP', denoise);
%                     re(re < 0) = 0;
%                     acc = sum(abs(re-s) <= tolerance) / n;
                    acc = snr(s,re);
                    res(loopcount, :) = [m, costratio, acc, noisel, sampling, snap];
                end
    %             acc = acc / snapn;

            end
            fprintf('sampling(%d)|m(%d)|k(%d) is over\n', sampling, m, k);
        end
        fprintf('sampling %d is over\n', sampling);
    end
end

for mi = 1 : mn
    m = ms(mi);
    res_ca(mi, :) = mean(res(res(:,1)==m & res(:,5)==0, 1:4));
    res_random(mi, :) = mean(res(res(:,1)==m & res(:,5)==1, 1:4));
    res_group(mi, :) = mean(res(res(:,1)==m & res(:,5)==2, 1:4));
    res_greedy(mi, :) = mean(res(res(:,1)==m & res(:,5)==3, 1:4));
end

figure
hold on
plot(res_random(:,2), res_random(:,3));
plot(res_group(:,2), res_group(:,3), 'r');
plot(res_ca(:,2), res_ca(:,3), 'm');
plot(res_greedy(:,2), res_greedy(:,3), 'k');

imethod = 'linear';
accstops = 0:0.01:1;
sres_target = [accstops' zeros(length(accstops), 1)];

for samplingi = 1 : samplingn
    sampling = samplings(samplingi);
    if (sampling == 0) tres = res_ca;
    elseif (sampling == 1) tres = res_random;
    elseif (sampling == 2) tres = res_group;
    elseif (sampling == 3) tres = res_greedy;
    end
    [sacc, idx] = sort(tres(:,3));
    [usacc, uidx] = unique(sacc);
    idx = idx(uidx);
    sres_target(:,2) = interp1([0 tres(idx,3)' 1], [0 tres(idx,2)' 1], accstops, imethod)';
    if (sampling == 0) sres_ca = sres_target;
    elseif (sampling == 1) sres_random = sres_target;
    elseif (sampling == 2) sres_group = sres_target;
    elseif (sampling == 3) sres_greedy = sres_target;
    end
end

% figure
% hold on
% plot(sres_random(:,1),sres_random(:,2))
% plot(sres_group(:,1),sres_group(:,2), 'r')
% plot(sres_ca(:,1),sres_ca(:,2), 'm')
% plot(sres_greedy(:,1),sres_greedy(:,2), 'k')

idx = [61 66 71 76 81 86 91 96];
% figure;
% bar([sres_random(idx,2) sres_group(idx,2) sres_ca(idx,2) sres_greedy(idx,2)]);
% sres_all = [sres_random(idx,2) sres_group(idx,2) sres_ca(idx,2) sres_greedy(idx,2)];
anoisel = mean(res(:,4));
[anoisel sres_random(idx,2)']
[anoisel sres_group(idx,2)']
% [anoisel sres_ca(idx,2)']
% [anoisel sres_greedy(idx,2)']