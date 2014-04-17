clear
clc

NDEBUG = 1;

workpath = 'e:\document\codes\matlab\wcs\';
week = 'w2';
seg = load('200seg.txt');
n = length(seg);
gt = load(['200gti_' week '.txt']);
% [intvn, ~] = size(gt);
% intvn = 26;
% positioncost = load('200positioncost_conn.txt');
positioncost = load('gpscost.txt');
% positioncost = ones(101,107);
dim = size(positioncost);
% timecost = load(['battercost_li_' week '.txt']);
timecost = load(['battercost_ex_' week '.txt']);
% timecost = load(['battercost_lo_' week '.txt']);
% timecost = ones(9000, intvn);
combine = 1; % 0 product, 1 sum;
recnum = load(['30Mrecnum_' week '.txt']);
% tN = sum(recnum(1:intvn));
cabcount = load(['cabcount_' week '.txt']);
cabID = find(cabcount  > 0);
cabnum = numel(cabID);

segmapping = zeros(dim(1), dim(2));
for pt = 1 : n
    segmapping(seg(pt,1) + 1, seg(pt,2) + 1) = pt;
end

recovery = 'adaptive-CS';

alpha_min = 0;
alpha_gap = 0.1;
alpha_max = 5;
alphas = alpha_min : alpha_gap : alpha_max;
alphas = [alphas 10 100];
alphan = numel(alphas);

mfrac_min = 0.30;
mfrac_gap = 0.05;
mfrac_max = 0.50;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracn = length(mfracs);

single_runs = mfracn * alphan;

intvns = [2,4,8,13,26,52,78,104,130];
% intvns = [130];
intvnn = numel(intvns);

runs = intvnn * single_runs;

%% load
all_rec = [];
stop_intv = max(intvns);
cache = 0;
for intv = 1 : stop_intv
    filename = sprintf('%s_intv%d.txt', week, intv);
    nrec = load([workpath 'taxidata\' filename]);
    all_rec = [all_rec; nrec];
    if (~NDEBUG)
        fprintf(' %d ', length(nrec));
        if (intv < stop_intv) fprintf('+'); end
    end
    if (~NDEBUG)
        if (intv - cache >= 10)
            fprintf('\n');
            cache = intv;
        end
    end
end
N = length(all_rec);

if (~NDEBUG)
    fprintf('= ');
end
fprintf('%d rec loaded\n', N);
cost = zeros(1, N);
sumcost = 0;
for i = 1 : N
    % generate cost
    if (combine == 0)
        cost(i) = positioncost(all_rec(i,4), all_rec(i,5)) * timecost(all_rec(i,1), all_rec(i,6));
    elseif (combine == 1)
        cost(i) = positioncost(all_rec(i,4), all_rec(i,5)) + timecost(all_rec(i,1), all_rec(i,6));
    end
    sumcost = sumcost + cost(i);
end
if (~NDEBUG) fprintf('cost loaded\n'); end

%%
exp_taxi_ca_adaptive;

%%
save('res_taxi_ca_adaptive.txt', 'res', '-ascii');
script_taxi_ca_adaptive_process_res;