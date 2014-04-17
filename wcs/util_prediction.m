workpath = 'e:\document\codes\matlab\wcs\';
week = 'w2';
seg = load('200seg.txt');
n = length(seg);
gt = load(['200gti_' week '.txt']);
[intvn, ~] = size(gt);
% intvn = 26;
% positioncost = load('200positioncost_conn.txt');
positioncost = load('gpscost.txt');
% positioncost = load('3gcost.txt');
%  positioncost = load('GPS+3G.txt');
% positioncost = ones(101,107);
dim = size(positioncost);
% timecost = load(['battercost_li_' week '.txt']);
% timecost = load(['battercost_ex_' week '.txt']);
% timecost = load(['battercost_lo_' week '.txt']);
timecost = ones(9000, intvn);
combine = 0;
recnum = load(['30Mrecnum_' week '.txt']);
tN = sum(recnum(1:intvn));
cabcount = load(['cabcount_' week '.txt']);
cabID = find(cabcount  > 0);
cabnum = numel(cabID);

segmapping = zeros(dim(1), dim(2));
for pt = 1 : n
    segmapping(seg(pt,1) + 1, seg(pt,2) + 1) = pt;
end
for i = 1 : 130
    file = sprintf('taxidata\\w2_intv%d.txt',i);
    nrec = load(file);
    fprintf('intv%d: %d recs\n', i, length(nrec));
end
cost = zeros(1, N);
sumcost = 0;
for i = 1 : N
    % generate cost
    if (combine == 0)
        cost(i) = positioncost(rec(i,4), rec(i,5)) * timecost(rec(i,1), rec(i,6));
    elseif (combine == 1)
        cost(i) = positioncost(rec(i,4), rec(i,5)) + timecost(rec(i,1), rec(i,6));
    end
    sumcost = sumcost + cost(i);
end

% -------------------------------------

alp=2.5;
thr = 0.859;
history_size = 10;
times = 1;
res_p = zeros(times, 1);
for timei = 1 : times
    cache = 0;
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
            wei = 1./(cintv*ones(1,historyn(cid)) - history_intv(cid, 1:historyn(cid)) + 1);
            pcost = history(cid, 1:historyn(cid))*wei'/sum(wei);
%             pcost = mean(history(cid, 1:historyn(cid)));
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
        end
        if (i - cache >= 1000000) 
            fprintf('%d processed\n', i);
            cache = i;
        end
    end
    res_p(timei) = mean(abs((cost-prediction)./cost));
end