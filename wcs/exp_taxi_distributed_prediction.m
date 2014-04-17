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
% 
% samplingn = length(samplings);
% mfracn = length(mfracs);
% runs = mfracn * samplingn;


%% to save

m_all = zeros(runs, intvn * n);
mask_all = m_all;
afrac = zeros(runs, 1);
sumcost = 0;
totalcost = afrac;
err = afrac;
acc = afrac;
ratiocost = afrac;
res = zeros(runs, 8);

intvprocessed = 0;
batchc = 0;
cache = 0;

while (intvprocessed < intvn)
    batchc = batchc + 1;
    
    if (~NDEBUG)
        fprintf('Batch work %d begin:\n', batchc);
    end
    
    rec = [];
    intv_l = (batchc - 1) * intvbatch + 1;
    intv_h = min(intvn, batchc * intvbatch);
    for intv = intv_l : intv_h
        filename = sprintf('%s_intv%d.txt', week, intv);
        nrec = load([workpath 'taxidata\' filename]);
        rec = [rec; nrec];
        if (~NDEBUG)
            fprintf(' %d ', length(nrec));
            if (intv < intv_h) fprintf('+'); end
        end
        if (~NDEBUG)
            if (intv - cache >= 10)
                fprintf('\n');
                cache = intv;
            end
        end
    end
    intvprocessed = intvprocessed + intvbatch;
    N = length(rec);
%     M = floor(N * mfrac);
    
    if (~NDEBUG)
        fprintf('= ');
    end
    fprintf('%d rec loaded\n', N);
    cost = zeros(1, N);
    acabcount = zeros(9000, 1);
    for i = 1 : N
        % generate cost
        if (combine == 0)
            cost(i) = positioncost(rec(i,4), rec(i,5)) * timecost(rec(i,1), rec(i,6));
        elseif (combine == 1)
            cost(i) = positioncost(rec(i,4), rec(i,5)) + timecost(rec(i,1), rec(i,6));
        end
        sumcost = sumcost + cost(i);
        acabcount(rec(i,1)) = acabcount(rec(i,1)) + 1;
    end
    if (~NDEBUG) fprintf('cost loaded\n'); end
    
    %% sampling
    run = 0;
    for samplingi = 1 : samplingn
        for mfraci = 1 : mfracn
            run = run + 1;
            sampling = samplings(samplingi);
            mfrac = mfracs(mfraci);
            res(run, 1) = sampling;
            M = floor(N * mfrac);
            eachM = floor(acabcount .* mfrac);
            samplex = zeros(1, N);
            
            tic;
            if (sampling == 0) % random sampling
                samplex = randperm(N, M);

            elseif (floor(sampling) == 1) % ca
                if (~exist('alpha', 'var') || alpha < 0) alpha = 1; end
                alp = alpha;
                thr = mfracs(mfraci);
                history_size = 10;
                history = zeros(cabnum, history_size);
                history_intv = zeros(cabnum, history_size);
                pointer = ones(cabnum, 1);
                historyn = zeros(cabnum, 1);
                prediction = zeros(1, N);
                for i = 1 : N
                    cid = rec(i, 1);
                    cintv = rec(i,6);
                    if (historyn(cid) <= 0)
                        key = 1;
                        pcost = cost(i);
                    else
%                         pcost = mean(history(cid, 1:historyn(cid)));
                        
                        wei = 1./(cintv*ones(1,historyn(cid)) - history_intv(cid, 1:historyn(cid)) + 1);
                        pcost = history(cid, 1:historyn(cid))*wei'/sum(wei);
                        pprob = pcost^(-alp)*50;
                        key = rand^(1/pprob);
                    end
                    prediction(i) = pcost;
                    if (key >= thr)
                        % pick
                        samplex(i) = 1;
                        history(cid, pointer(cid)) = cost(i);
                        history_intv(cid, pointer(cid)) = rec(i, 6);
                        pointer(cid) = pointer(cid) + 1;
                        if (pointer(cid) > history_size) pointer(cid) = 1; end
                        if (historyn(cid) < history_size) historyn(cid) = historyn(cid) + 1; end
%                         fprintf('rec %d picked: key=%f cost=%f prediction=%f\n', i, key, cost(i), pcost);
                        
                    else
%                         fprintf('rec %d not picked: key=%f cost=%f prediction=%f\n', i, key, cost(i), pcost);
                    end
                end
                samplex = find(samplex==1);

            elseif (floor(sampling) == 2) %pair
                group = round(10 * (sampling - floor(sampling)));
                if (group == 0) group = 2; end
                gs = randperm(N);
                gn = min(floor(N/group), M);
                sn = floor(M/gn) * ones(1, gn);
                toadd = M - gn*floor(M/gn);
                if (toadd > 0)
                    sn(1:toadd) = sn(1:toadd) + 1;
                end
                chosen = 0;
                for gi = 1 : gn
                    [~, midx] = sort(cost(gs((gi-1)*group+1 : gi*group)));
                    samplex(chosen+1:chosen+sn(gi)) = gs((gi-1)*group*ones(1, sn(gi))+midx(1:sn(gi)));
                    chosen = chosen + sn(gi);
                end
%                 p2 = randperm (N, group * M);
%                 for chosen = 1 : M
%                     [~, midx] = min(cost(p2((chosen-1)*group+1 : chosen*group)));
%                     samplex(chosen) = p2((chosen-1)*group + midx);
%                 end

            elseif (sampling == 3) % ia
                u = 10;
                ipc = ones(1, N) ./ cost;
                for chosen = 1 : M
                    cands = find(ipc==max(ipc));
                    samplex(chosen) = cands(randsample(length(cands), 1));
                    % fading
                    chosen_rec = rec(samplex(chosen), :);
                    xvec = (chosen_rec(2) .* ones(1,N) - rec(:,2)').^2;
                    yvec = (chosen_rec(3) .* ones(1,N) - rec(:,3)').^2;
                    distvec = -1 .* (xvec + yvec).^(-u/2);
                    fade = u.^distvec;
                    fade(samplex(chosen)) = 0;
                    ipc = ipc .* fade;
                end
            elseif (sampling == 4) % TMC
                cabM = floor(mfrac * cabnum);
                chosenCab = zeros(1, max(cabID));
                chosenCab(cabID(randsample(cabnum, cabM))) = ones(1, cabM);
                samplex = [];
                samc = 0;
                for i = 1 : N
                    if (chosenCab(rec(i, 1)) == 1)
                        samc = samc + 1;
                        samplex(samc) = i;
                    end
                end
            end
            res(run, 8) = toc;

            % remove duplicate
            samplex = unique(samplex);
            aM = numel(samplex);
            counter = zeros(intvn, n);
%             speeds = -1*ones(intvn*n, 2);
            for i = 1 : aM
                segidx = segmapping(rec(samplex(i),4) + 1, rec(samplex(i),5) + 1);
                if (segidx > 0)
                    intv = rec(samplex(i), 6);
                    totalcost(run) = totalcost(run) + cost(samplex(i));
                    counter(intv, segidx) = counter(intv, segidx) + 1;
                    m_all(run, (segidx-1)*intvn + intv) = m_all(run, (segidx-1)*intvn + intv) + rec(samplex(i), 7);
%                     speeds((segidx-1)*intvn + intv, counter(intv, segidx)) = rec(samplex(i), 7);
                end
            end
            afrac(run) = afrac(run) + sum(sum(counter));
            for intv = intv_l : intv_h
                for i = 1 : n
                    if (counter(intv, i) > 1)
                        m_all(run, (i-1)*intvn + intv) = m_all(run, (i-1)*intvn + intv) / counter(intv, i);
%                         speed_ok = sort(speeds((i-1)*intvn + intv, 1:counter(intv,i)));
%                         cut = floor(0.1*counter(intv,i));
%                         if (cut == 0) cut = 1; end
%                         m_all(run, (i-1)*intvn + intv) = mean(speed_ok(cut + 1:(counter(intv,i)-cut)));
                        mask_all(run, (i-1)*intvn + intv) = 1;
                    else m_all(run, (i-1)*intvn + intv) = 0; 
                    end
                end
            end
            fprintf('Batch work %d sampling(%.1f) finished: %d recs(%.2f) sampled\n', batchc, sampling, sum(sum(counter)), mfrac);
        end
    end
    fprintf('Batch work %d finished\n', batchc);
end
afrac = afrac ./ tN;
coverage = sum(mask_all, 2) ./ intvn ./ n;

%% recovery
run = 0;
for samplingi = 1 : samplingn
    for mfraci = 1 : mfracn
        run = run + 1;
        mask = zeros(intvn, n);
        measure = mask;
        mask(:) = mask_all(run,:);
        measure(:) = m_all(run,:);
        if (strcmp(recovery, 'SVT'))
           %% SVT
            n1 = intvn; n2 = n; rk = 5;    
            df = rk*(n1+n2-rk);
            oversampling = 3; 
            p  = afrac(run);

            fprintf('Matrix completion: %d x %d matrix, rank %d, %.1f%% observations\n',...
                n1,n2,rk,100*p);
            tau = oversampling*sqrt(n1*n2);
            delta = 1.2/p;   
            %{
             if n1 and n2 are very different, then
               tau should probably be bigger than 5*sqrt(n1*n2)

             increase tau to increase accuracy; decrease it for speed

             if the algorithm doesn't work well, try changing tau and delta
               i.e. if it diverges, try a smaller delta (e.g. delta < 2 is a 
               safe choice, but the algorithm may be slower than necessary).
            %}
            maxiter = 500; 
            tol = 1e-4;
            % Note: SVT, as called below, is setup for noiseless data 
            %   (i.e. equality constraints).

            fprintf('\nSolving by SVT...\n');
            tic
            idx  = find(mask==1);
            [U,S,V,numiter] = SVT([n1 n2],idx,measure(idx),tau,delta,maxiter,tol);
            %         [U,S,V,numiter] = FPC([n1 n2],idx,gt(idx),0.01,maxiter,0.1);
            toc

            r = U*S*V';

            % Show results
            fprintf('The recovered rank is %d\n',length(diag(S)) );
            fprintf('The relative error on Omega is: %d\n', norm(gt(idx)-r(idx))/norm(gt(idx)))
            fprintf('The relative recovery error is: %d\n', norm(gt-r,'fro')/norm(gt,'fro'))
            fprintf('The relative recovery in the spectral norm is: %d\n', norm(gt-r)/norm(gt))

            elseif (strcmp(recovery, 'INEXACT'))
           %%
                entries = zeros(sum(sum(mask)), 3);
                enc = 0;
                for i = 1 : intv
                    for j = 1 : n
                        if (mask(i,j)==1)
                            enc = enc + 1;
                            entries(enc,1:3) = [i,j,measure(i,j)];
                        end
                    end
                end
                m_spm = spconvert(entries);
                [R, iter, rk] = inexact_alm_mc(m_spm, 0.1);
                r = R.U * R.V';
                r(find(r < 0)) = zeros(1, numel(find(r < 0)));

            elseif (strcmp(recovery, 'MC'))
                %% MC
                rk = 20;
                lam = 32;
                r = matrixcompletion(measure, mask, rk, lam, 300);
                r(find(r < 0)) = zeros(1, numel(find(r < 0)));

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
                        fprintf('sampling(%.1f | %.2f) recovery finished: k = %d, round %d\n',samplings(samplingi), mfracs(mfraci), k, i);
                    end
                end
                r(find(r < 0)) = zeros(1, numel(find(r < 0)));
        
        else
        %% CS snapshot
        %     r = zeros(intvn, n);
        %     [b,~] = princomp(gt);
        %     b = b';
        %     k = floor(m/log(n));
        %     for intv = 1 : intvn
        %         samplex = find(mask_all(intv, :)==1);
        %         y = m_all(intv, samplex);
        %         r(intv, :) = my_csf(y', samplex, b, k, 'OMP', 1)';
        %         r(find(r < 0)) = zeros(1, numel(find(r < 0)));
        %         %error = sum(sum(abs((gt(intv, :)-r(intv, :)).*unavail(intv, :))))/sum(sum(abs(gt(intv, :).*unavail(intv, :))));
        %         %fprintf('recovery finished intv %d: %f\n', intv, error);
        %     end
        end
        
       %% res
        erm = abs(gt(1:intvn,1:n)-r);
        err(run) = sum(sum(abs(gt(1:intvn,1:n)-r)))/sum(sum(abs(gt(1:intvn,1:n))));
        acc(run) = numel(find(erm <= 10))/intv/n;
        ratiocost(run) = totalcost(run) / sumcost;
        fprintf('sampling=%.1f frac=%f acc=%f err=%f coverage=%f | cost=%f ratio=%f\n', res(run,1), afrac(run), acc(run), err(run), coverage(run), totalcost(run), ratiocost(run));
    end
end
res(:, 2:7) = [afrac, acc, err, coverage, totalcost, ratiocost];