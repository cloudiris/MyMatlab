%---- finding optimal m-samples in 2D scenario ------

clear
clc

%% parameter setup

workpath = 'd:\document\codes\matlab\wcs\';
datapath = 'd:\document\codes\data\';
data = 'temp.grid.64.64';
s2 = load([datapath data '.txt']);
s = hilbertcurve(s2)';

si = size(s2);
n = si(1) * si(2);
k = floor(0.01 * n);
b = base_dct4(n);
method = 'OMP';

mmin = 410;
mgap = 410;
mmax = 1230;
looptime = floor((mmax - mmin) / mgap) + 1;
mstr = sprintf('%d-%d-%d', mmin, mgap, mmax);

%%% cost distribution
% cost_distribution = 'fixed';
% cost = load([datapath 'cost\cost.' cost_distribution '.txt']);
% cost = generatecost(32, 32, cost_distribution);
cost = load('exp1.txt');
costv = hilbertcurve(cost);
allcost = sum(costv);

%%% measurement
acc_measurement = 'errpluscost';
if (strcmp(acc_measurement, 'snr'))
    accfunc = @snr;
    argfunc = @divide;
    comp = @greater;
    accidx = 2;
elseif (strcmp(acc_measurement, 'mse'))
    accfunc = @mse;
    comp = @less;
    accidx = 3;
elseif (strcmp(acc_measurement, 'acc100'))
    accfunc = @acc100;
    comp = @greater;
    accidx = 4;
elseif (strcmp(acc_measurement, 'l2norm'))
    accfunc = @l2normratio;
    comp = @greater;
    accidx = 2;
elseif (strcmp(acc_measurement, 'errpluscost'))
    accfunc = @err01;
    argfunc = @errpluscost;
    comp = @less;
    accidx = 2;
end

%%% save filename
% savefile = ['optimal(' acc_measurement ').data(' data ').cost(' cost_distribution ').mat'];
% savesamplex = ['optimal(samplex).data(' data ').cost(' cost_distribution ').mat' ];
savefile = '2d_optimal.txt';
savesamplex = '2d_optimal_samplex.txt';
target = 'globaloptimal';

if (exist([workpath savefile], 'file'))
    load([workpath savefile]);
else
    globaloptimal = inf * ones(n, 7);
end

if (exist([workpath savesamplex], 'file'))
    load([workpath savesamplex]);
else
    optimalsamplex = zeros(n);
end

%% global annealing parameter
iteration = 100000;

%% main loop
for mcount = 1 : looptime
    %% initialization
    m = mmin + (mcount - 1) * mgap;
    
%     randp = randperm(n);
%     samplex = randp(1 : m);
    [~, in] = sort(costv);
    samplex = in(1 : m);
    % greedy

    samples = s(samplex);
    totalcost = sum(costv(samplex));
    costratio = totalcost / allcost;
    re = my_csf(samples, samplex, b, k, method, 1);
%     cacc = accfunc(s, re);
    cacc = numel(find(abs((s-re)./s)<=0.01)) / numel(s);
    optimal = cacc/costratio;
    if (globaloptimal(m, 1) == inf || comp(optimal, globaloptimal(m, 1)))
        globaloptimal(m, 1) = optimal;
        globaloptimal(m, 2) = snr(s, re);
        globaloptimal(m, 3) = mse(s, re);
        globaloptimal(m, 4) = numel(find(abs((s-re)./s)<=0.01)) / numel(s);
        globaloptimal(m, 5) = norm(re, 2);
        globaloptimal(m, 6) = totalcost;
        globaloptimal(m, 7) = costratio;
        optimalsamplex(m, :) = zeros(1, n);
        optimalsamplex(m, samplex) = ones(1, m);
    end
    
    output = sprintf('m = %d: Iteration ', m);
    loopclearline = numel(output);
    fprintf(output);
    output = sprintf('%d/%d', mcount, iteration);
    fprintf(output);
    clearline = numel(output);

    %% annealing
    %%% annealing parameters setup
%     iteration = 10000;
    temperature = 500;
    cooling = 0.99;
    LMD = 1; 

    %%% records
%     changes = zeros(iteration, 2);
%     originalSamplex = samplex;
%     res = zeros(iteration, 3);
%     res(1, :) = [cacc totalcost optimal];

    % i = 1;
    for i = 2 : iteration
        %% cooling
        temperature = temperature * cooling;

        %% random transition
        %%% randomly revome one sample from *samplex*
        %%% randomly add one sample to *samplex*
        randp = randperm(m);
        removeindex = randp(1);
        removed = samplex(removeindex);
        one2n = 1 : n;
        one2n(samplex) = 0;
        randp = randperm(n);
        j = 1;
        while (one2n(randp(j)) == 0)
            j = j + 1;
        end
        added = randp(j);
        samplex(removeindex) = added;

        %% update cost
        samples = s(samplex);
        totalcost = sum(costv(samplex));
        costratio = totalcost / allcost;
        re = my_csf(samples, samplex, b, k, method, 1);
        cacc = numel(find(abs((s-re)./s)<=0.01)) / numel(s);
        cr = cacc/costratio;
        if (comp(cr, optimal))
            %% a better case, ACCEPT
%             changes(i, 1) = removed;
%             changes(i, 2) = added;
%             res(i, 1) = cacc;
%             res(i, 2) = costratio;
%             res(i, 3) = cr;
            optimal = cr;
            if (comp(optimal, globaloptimal(m, 1)))
                globaloptimal(m, 1) = optimal;
                globaloptimal(m, 2) = snr(s, re);
                globaloptimal(m, 3) = mse(s, re);
                globaloptimal(m, 4) = numel(find(abs((s-re)./s)<=0.01)) / numel(s);
                globaloptimal(m, 5) = norm(re, 2);
                globaloptimal(m, 6) = totalcost;
                globaloptimal(m, 7) = costratio;
                optimalsamplex(m, :) = zeros(1, n);
                optimalsamplex(m, samplex) = ones(1, m);
            end
            %% print info
%             clc;
            fprintf(repmat('\b', 1, clearline));
            output = sprintf('%d/%d', i, iteration);
            fprintf(output);
            clearline = numel(output);
    
        else
            %% a no-better case, accept w/ probability
            pr = exp(cr - optimal) / LMD / temperature;
            acc = rand();
            if (acc <= pr)
                optimal = cr;
            else
                %% reject
                %%% restore transition
                samplex(removeindex) = removed;
                temperature = temperature / cooling;
            end
            %% print info
%             clc;
            fprintf(repmat('\b', 1, clearline));
            output = sprintf('%d/%d', i, iteration);
            fprintf(output);
            clearline = numel(output);
            
        end
        %% iteration figure update
%         figure(1);
%         hold on;
%         plot(i, optimal);
    end
    fprintf(repmat('\b', 1, clearline));
    fprintf(repmat('\b', 1, loopclearline));
    fprintf('m = %d : Acc = %.3f, Cost = %.3f%%, Optimal = %.3f\n', m, globaloptimal(m, accidx), globaloptimal(m, 7)*100, globaloptimal(m, 1));
%     clf;
%     figure(1);
%     plot(globaloptimal(1:mcount, 1), globaloptimal(1:mcount, 2));
end

save(savefile, target);
save(savesamplex, 'optimalsamplex');