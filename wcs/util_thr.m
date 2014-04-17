% b = load('battercost_ex_w2.txt');
clear;
% gpscost = load('gpscost.txt');
% gpscost = load('gps+3g.txt');
gpscost = load('3gcost.txt');
rec = [];
for i = 1 : 130
    file = sprintf('taxidata\\w2_intv%d.txt',i);
    r = load(file);
    rec = [rec; r];
    fprintf('%d\n',i);
end
N = length(rec);
cost = zeros(1, N);
for i = 1 : N
    lat = rec(i,4);
    lon = rec(i,5);
    cid = rec(i,1);
    intv = rec(i,6);
    cost(i) = gpscost(lat, lon);
end

alps = [2.5];
fac = 50;
for i = 1 : numel(alps)
    alp = alps(i);
    w = cost.^(-alp)*fac;
    key = (rand(1, N)).^(1./w);
    thr(i,:) = prctile(key,80:-10:50);
end
fprintf('%f %f %f %f', thr(1,1), thr(1,2), thr(1,3), thr(1,4));
% save('thr.txt','thr','-ascii');