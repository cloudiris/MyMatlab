% clear
% clc

NDEBUG = 1;

workpath = 'e:\document\codes\matlab\wcs\';
week = 'w2';
seg = load('200seg.txt');
n = length(seg);
gt = load(['200gti_' week '.txt']);
[intvn, ~] = size(gt);
% intvn = 26;
% positioncost = load('200positioncost_conn.txt');
% positioncost = load('gpscost.txt');
% positioncost = load('3gcost.txt');
positioncost = load('gps+3g.txt');
% positioncost = ones(101,107);
dim = size(positioncost);
% timecost = load(['battercost_li_' week '.txt']);
% timecost = load(['battercost_ex_' week '.txt']);
% timecost = load(['battercost_lo_' week '.txt']);
timecost = ones(9000, intvn);
combine = 0; % 0 product, 1 sum;
recnum = load(['30Mrecnum_' week '.txt']);
tN = sum(recnum(1:intvn));
cabcount = load(['cabcount_' week '.txt']);
cabID = find(cabcount  > 0);
cabnum = numel(cabID);

segmapping = zeros(dim(1), dim(2));
for pt = 1 : n
    segmapping(seg(pt,1) + 1, seg(pt,2) + 1) = pt;
end

recovery = 'CS';
intvbatch = intvn;

ks = [1,2,3,4,5,10,20];
kn = numel(ks);

mfrac_min = 0.05;
mfrac_gap = 0.05;
mfrac_max = 0.9;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracn = length(mfracs);

times = 1;
single_runs = mfracn * times;

runs = kn * single_runs;

exp_taxi_k;

save('res_taxi_k.txt', 'res', '-ascii');

tointerp = 0.7 : 0.05 : 0.95;
res_k = zeros(kn + 1, numel(tointerp));
res_k(kn + 1, :) = tointerp;
f = fopen('res_k_fixacc_cost.txt', 'w');
fprintf(f, 'k ');
for i = 1 : numel(tointerp)
    fprintf(f, '%f ', tointerp(i));
end
fprintf(f, '\n');
for ki = 1 : kn
    k = ks(ki);
    res_tmp = res(find(res(:,1) == k), 2:8);
    x = res_tmp(:,2);
    y = res_tmp(:,6);
    icost = interp1(x, y, tointerp, 'spline');
    res_k(ki, :) = icost;
    fprintf(f, '%f ', k);
    for i = 1 : numel(icost)
        fprintf(f, '%f ', icost(i));
    end
    fprintf(f, '\n');
end
fclose(f);
plot(ks, res_k(1:kn, :)');
for i = 1 : numel(tointerp)
    l(i, :) = sprintf('acc=%1.2f', tointerp(i));
end
legend(l);
