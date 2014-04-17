% clear
% clc
% NDEBUG = 0;
%% pretreatment
% N = 10000;
% low = 100;
% high = 200;
% 
% cost = low + (high - low) .* rand(1, N);
% % cost = normrnd(3,1, 1,N);
% 
% cost = cost(find(cost > low));
% cost = cost(find(cost < high));
% N = numel(cost);
% sumcost = sum(cost);

%% parameter
% 0 random; 1 ca; 2 pair; 3 ia; 4 TMC
% if (~exist('samplings', 'var')) samplings = 0; end
% if (~exist('mfracs', 'var')) mfracs = 0.01; end
% recovery = 'CS';
% 
% samplingn = length(samplings);
% mfracn = length(mfracs);
% runs = mfracn * samplingn;


%% to save

mask_all = zeros(runs, N);
afrac = zeros(runs, 1);
totalcost = afrac;
ratiocost = afrac;
res = zeros(runs, 6);

cache = 0;

%% sampling
run = 0;
for samplingi = 1 : samplingn
    for mfraci = 1 : mfracn
        for timei = 1 : times
            run = run + 1;
            sampling = samplings(samplingi);
            mfrac = mfracs(mfraci);
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
                res(run, 2) = alp;
                prob = cost.^(-1 * alp) / sum(cost.^(-1 * alp));
                samplex = randsamplewtr(N, M, prob);

            elseif (floor(sampling) == 2) %pair
                group = round(10 * (sampling - floor(sampling)));
                if (group == 0) group = 2; end
                res(run, 2) = group;
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
    %         elseif (sampling == 4) % TMC
    %             cabM = floor(mfrac * cabnum);
    %             chosenCab = zeros(1, max(cabID));
    %             chosenCab(cabID(randsample(cabnum, cabM))) = ones(1, cabM);
    %             samplex = [];
    %             samc = 0;
    %             for i = 1 : N
    %                 if (chosenCab(rec(i, 1)) == 1)
    %                     samc = samc + 1;
    %                     samplex(samc) = i;
    %                 end
    %             end
            end
            res(run, 6) = toc;
            mask_all(run, samplex) = 1;

            % remove duplicate
            samplex = unique(samplex);
            aM = numel(samplex);
            afrac(run) = aM;
            totalcost(run) = sum(cost(samplex));
            ratiocost(run) = totalcost(run) / sumcost;
            fprintf('sampling(%.1f) finished: frac(%.2f) | cost(%.5f)\n', sampling, mfrac, ratiocost(run));
        end
    end
end

afrac = afrac ./ N;
res(:, 3:5) = [afrac, totalcost, ratiocost];