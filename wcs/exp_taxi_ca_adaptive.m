% clear
% clc
% NDEBUG = 0;
%% pretreatment
% workpath = 'd:\document\codes\matlab\wcs\';
% week = 'w2';
% seg = load('200seg.txt');
% n = length(seg);
% gt = load(['200gti_' week '.txt']);
% [intvn, ~] = size(gt);
% positioncost = load('200positioncost_station.txt');
% % positioncost = load('test_positioncost.txt');
% % positioncost = ones(101,107);
% % positioncost = load('200positioncost_uniform.txt');
% dim = size(positioncost);
% timecost = ones(1, intvn);
% % timecost = load(['30Mtimecost_' week '.txt']);
% % timecost = timecost + noise_n(intvn, 0.1*norm(timecost))'; % add noise
% recnum = load(['30Mrecnum_' week '.txt']);
% tN = sum(recnum);
% cabcount = load(['cabcount_' week '.txt']);
% cabID = find(cabcount  > 0);
% cabnum = numel(cabID);
% 
% segmapping = zeros(dim(1), dim(2));
% for pt = 1 : n
%     segmapping(seg(pt,1) + 1, seg(pt,2) + 1) = pt;
% end

%% parameter
% 0 random; 1 ca; 2 pair; 3 ia; 4 TMC
% if (~exist('samplings', 'var')) samplings = 0; end
% if (~exist('mfracs', 'var')) mfracs = 0.01; end
% intvbatch = intvn;
% recovery = 'CS';


%% to save

res = zeros(runs, 9);
afrac = zeros(runs, 1);
totalcost = afrac;
err = afrac;
acc = afrac;
ratiocost = afrac;
coverage = afrac;

run = 0;
for intvni = 1 : intvnn
    intvn = intvns(intvni);
    rec = all_rec(find(all_rec(:,6) <= intvn),:);
    N = length(rec);
    sumcost = sum(cost(1:N));    

%% sampling
    for alphai = 1 : alphan
        alpha = alphas(alphai);

        for mfraci = 1 : mfracn
            mfrac = mfracs(mfraci);

            run = run + 1;
            res(run, 1) = alpha;
            res(run, 9) = intvn;
            M = floor(N * mfrac);
            samplex = zeros(1, M);

            tic;

            if (alpha < 0) alpha = 0; end
            alp = alpha;
            prob = cost.^(-1 * alp) / sum(cost.^(-1 * alp));
            samplex = randsamplewtr(N, M, prob(1:N));

            res(run, 8) = toc;
            
            % remove duplicate
            samplex = unique(samplex);
            aM = numel(samplex);
            counter = zeros(intvn, n);
            m_all = zeros(intvn, n);
            mask_all = m_all;
            for i = 1 : aM
                segidx = segmapping(rec(samplex(i),4) + 1, rec(samplex(i),5) + 1);
                if (segidx > 0)
                    intv = rec(samplex(i), 6);
                    totalcost(run) = totalcost(run) + cost(samplex(i));
                    counter(intv, segidx) = counter(intv, segidx) + 1;
                    m_all(intv, segidx) = m_all(intv, segidx) + rec(samplex(i), 7);
                end
            end
            afrac(run) = afrac(run) + sum(sum(counter));
            afrac(run) = afrac(run) / N;
            for intv = 1 : intvn
                for i = 1 : n
                    if (counter(intv, i) > 1)
                        m_all(intv, i) = m_all(intv, i) / counter(intv, i);
    %                         speed_ok = sort(speeds((i-1)*intvn + intv, 1:counter(intv,i)));
    %                         cut = floor(0.1*counter(intv,i));
    %                         if (cut == 0) cut = 1; end
    %                         m_all(run, (i-1)*intvn + intv) = mean(speed_ok(cut + 1:(counter(intv,i)-cut)));
                        mask_all(intv, i) = 1;
                    else m_all(intv, i) = 0; 
                    end
                end
            end
            fprintf('sampling(%.1f) finished: %d recs(%.2f) sampled\n', alpha, sum(sum(counter)), mfrac);
            
            %% recovery
            mask = mask_all;
            measure = m_all;
            coverage(run) = sum(sum(mask)) / intvn / n;
            if (strcmp(recovery, 'SVT'))
              
            elseif (strcmp(recovery, 'INEXACT'))
              
            elseif (strcmp(recovery, 'MC'))
               
            elseif (strcmp(recovery, 'CS'))
           %% CS time-series
                r = zeros(intvn, n);
                nn = 10 * intvn;
                base = base_dct4(nn);
                for i = 1 : 211
                    samplex = find(mask((i - 1) * nn + 1 : i * nn)==1);
                    if (isempty(samplex)) continue; end
                    y = measure(samplex + (i-1)*nn);
                    k = round(length(samplex) / log2(nn) / 2.5);
                    if (k < 1) k = 1; end
                    r((i - 1) * nn + 1 : i * nn) = my_csf(y', samplex, base, k, 'OMP', 1);

                    if (~NDEBUG)
                        fprintf('alpha=(%.1f | %.2f) recovery finished: k = %d, round %d\n',alphas(alphai), mfracs(mfraci), k, i);
                    end
                end
                r(find(r < 0)) = zeros(1, numel(find(r < 0)));

            elseif (strcmp(recovery, 'adaptive-CS'))
                r = zeros(intvn, n);
                nn = 1300;
                col = ceil(nn / intvn);
                nn = col * intvn;
                rounds = ceil(n / col);
                eachround = nn * ones(1, rounds);
                eachround(rounds) = n * intvn - sum(eachround(1:rounds-1));
                base = base_dct4(nn);
                for i = 1 : rounds
                    if (i == rounds) base = base_dct4(eachround(i)); end
                    start = sum(eachround(1:i-1));
                    range = start + 1 : start + eachround(i);
                    samplex = find(mask(range)==1);
                    if (isempty(samplex)) continue; end
                    y = measure(samplex + start);
                    k = round(length(samplex) / log2(nn) / 2.5);
                    if (k < 1) k = 1; end
                    r(start + 1 : start + eachround(i)) = my_csf(y', samplex, base, k, 'OMP', 1);

                    if (~NDEBUG)
                        fprintf('alpha=(%.1f | %.2f) recovery finished: k = %d, round %d\n',alphas(alphai), mfracs(mfraci), k, i);
                    end
                end
                r(find(r < 0)) = zeros(1, numel(find(r < 0)));
            else
          
            end

           %% res
            erm = abs(gt(1:intvn,1:n)-r);
            err(run) = sum(sum(abs(gt(1:intvn,1:n)-r)))/sum(sum(abs(gt(1:intvn,1:n))));
            acc(run) = numel(find(erm <= 10))/intv/n;
            ratiocost(run) = totalcost(run) / sumcost;
            fprintf('sampling=%.1f frac=%f acc=%f err=%f coverage=%f | cost=%f ratio=%f\n', res(run,1), afrac(run), acc(run), err(run), coverage(run), totalcost(run), ratiocost(run));
        end
    end
%     coverage = sum(mask_all, 2) ./ intvn ./ n;
end

res(:, 2:7) = [afrac, acc, err, coverage, totalcost, ratiocost];