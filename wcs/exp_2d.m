% clear
% clc
% NDEBUG = 0;
%% pretreatment
% workpath = 'd:\document\codes\matlab\wcs\';
% datapath = 'd:\document\codes\data\';
% data = 'temp.grid.64.64.txt';
% s2 = load([datapath data]);
% s = hilbertcurve(s2)';
% N = numel(s);
% 
% positioncost = load('norm1.txt');
% % positioncost = load('norm2.txt');
% % positioncost = load('norm3.txt');
% 
% cost = hilbertcurve(positioncost);
% sumcost = sum(cost);
% b = base_dct4(N);

%% parameter
% 0 random; 1 ca; 2 pair; 3 ia; 4 TMC
if (~exist('samplings', 'var')) samplings = 1; end
if (~exist('mfracs', 'var')) mfracs = 0.01; end
if (~exist('times', 'var')) times = 10; end
recovery = 'CS';

samplingn = length(samplings);
mfracn = length(mfracs);
single_runs = mfracn * times;
runs = single_runs * samplingn;


%% to save

mask_all = zeros(runs, N);
afrac = zeros(runs, 1);
totalcost = afrac;
err = afrac;
acc = afrac;
coverage = afrac;
ratiocost = afrac;
res = zeros(runs, 8);

%% sampling
run = 0;
for samplingi = 1 : samplingn
    sampling = samplings(samplingi);
    
    for mfraci = 1 : mfracn
        mfrac = mfracs(mfraci);
        
        for timei = 1 : times
            run = run + 1;  
            res(run, 1) = sampling;
            M = floor(N * mfrac);
            samplex = zeros(1, M);

            tic;
            if (sampling == 0) % random sampling
                samplex = randperm(N, M);

            elseif (floor(sampling) == 1) % ca
                if (sampling == 1.1) alp = 0.1;
                elseif (sampling == 1.2) alp = 0.5;
                elseif (sampling == 1.3) alp = 2;
                elseif (sampling == 1.4) alp = 5;
                elseif (sampling == 1.5) alp = 100;
                else alp = 1; end
                prob = cost.^(-1 * alp);
                samplex = randsamplewtr(N, M, prob);

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
                
            elseif (sampling == 5) % greedy
                [~, cheap] = sort(cost);
                samplex = cheap(1:M);
                
            elseif (sampling == 6) % local group
                radius = 5;
                mask = zeros(dim);
                bkc = positioncost;
                for i = 1 : M
                    start = [randi([1 (dim(1) - radius + 1)]), randi([1 (dim(2) - radius + 1)])];
                    while (numel(find(bkc(start(1):start(1)+radius-1, start(2):start(2)+radius-1)<inf)) == 0) 
                        start = [randi([1 (dim(1) - radius + 1)]), randi([1 (dim(2) - radius + 1)])]; 
                    end
                    [mvs, ins] = min(bkc(start(1):start(1) + radius - 1, start(2):start(2)+radius - 1));
                    [~, y] = min(mvs);
                    x = ins(y);
                    x = start(1) + x-1;
                    y = start(2) + y-1;
                    mask(x,y) = 1;
                    bkc(x,y) = inf;
                end
                samplex = zeros(1, M);
                maskv = hilbertcurve(mask);
                samplex(:) = find(maskv(:)==1);
                
            end
            res(run, 8) = toc;

            % remove duplicate
            samplex = unique(samplex);
            totalcost(run) = sum(cost(samplex));
            ratiocost(run) = totalcost(run) / sumcost;
            aM = numel(samplex);
            afrac(run) = aM/N;
            
            % recovery
            k = round(aM / log2(N) / 2.5);
            if (k < 1) k = 1; end
            r = my_csf(s(samplex), samplex, b, k, 'OMP', 1);
            
            e = s - r;
            acc(run) = snr(s, r);
            err(run) = norm(s-r, 1) / norm(s, 1);
            coverage(run) = numel(find(abs(e(2049:4096)./s(2049:4096)) <= 0.03))/N*2;
            fprintf('Run %d sampling(%.1f) finished: %.3f sampled | correctness = %.3f | cost= %.3f\n', run, sampling, afrac(run), coverage(run), ratiocost(run));
        end
    end
end

res(:, 2 : 7) = [afrac, acc, err, coverage, totalcost, ratiocost];