clear
clc

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
% positioncost = load('3Gcost.txt');
positioncost = load('GPS+3G.txt');
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

alpha_min = 0;
alpha_gap = 0.5;
alpha_max = 5;
alphas = alpha_min : alpha_gap : alpha_max;
alphas = [alphas 10 100];
% alphas = [5 10 100];
alphan = numel(alphas);

mfrac_min = 0.05;
mfrac_gap = 0.05;
mfrac_max = 0.95;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracn = length(mfracs);

times = 1;
single_runs = mfracn * times;

runs = alphan * single_runs;

exp_taxi_alpha;

save('res_taxi_alpha.txt', 'res', '-ascii');

tointerp = 0.7 : 0.05 : 0.95;
res_alpha = zeros(alphan + 1, numel(tointerp));
res_alpha(alphan + 1, :) = tointerp;
f = fopen('res_alpha_fixacc_cost.txt', 'w');
fprintf(f, 'alpha ');
for i = 1 : numel(tointerp)
    fprintf(f, '%f ', tointerp(i));
end
fprintf(f, '\n');
for alphai = 1 : alphan
    alpha = alphas(alphai);
    res_tmp = res(find(res(:,1) == alpha), 2:8);
    x = res_tmp(:,2);
    y = res_tmp(:,6);
    icost = interp1(x, y, tointerp, 'spline');
    res_alpha(alphai, :) = icost;
    fprintf(f, '%f ', alpha);
    for i = 1 : numel(icost)
        fprintf(f, '%f ', icost(i));
    end
    fprintf(f, '\n');
end
fclose(f);
plot(alphas, res_alpha(1:alphan, :)');
grid on;
set(gca,'linewidth',2);
for i = 1 : numel(tointerp)
    l(i,:) = sprintf('acc=%1.2f', tointerp(i));
end
legend(l);