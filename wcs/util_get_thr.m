
for i = 1 : length(rates)
    intvn = rates(i,1);
    rate = rates(i,2);
    alp = rates(i,3);
    partialcost = cost(1:sum(recnum(1:intvn)));
    partialN = numel(partialcost);
    w = partialcost.^(-alp)*50;
    key = rand(1,partialN).^(1./w);
    thr = prctile(key, (1-rate)*100);
    rates(i, 4) = thr;
    fprintf('intvn=%d, rate=%f, alp=%f, thr=%f\n', intvn, rate, alp, thr);
end