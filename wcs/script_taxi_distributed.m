clear
clc

NDEBUG = 0;

workpath = 'd:\document\codes\matlab\wcs\';
week = 'w2';
seg = load('200seg.txt');
n = length(seg);
gt = load(['200gti_' week '.txt']);
[intvn, ~] = size(gt);
% intvn = 26;
% positioncost = load('200positioncost_conn.txt');
positioncost = load('gpscost.txt');
% positioncost = ones(101,107);
dim = size(positioncost);
% timecost = load(['battercost_li_' week '.txt']);
% timecost = load(['battercost_ex_' week '.txt']);
% timecost = load(['battercost_lo_' week '.txt']);
% timecost = ones(9000, intvn);
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
mfrac_min = 0.5;
mfrac_gap = 0.05;
mfrac_max = 0.5;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracs = [ 0.9951    0.9847    0.9654];
mfracn = length(mfracs);

times = 1;
single_runs = mfracn * times;

samplings = [1];
samplingn = length(samplings);

cest = 3e8;

runs = samplingn * mfracn;

exp_taxi_distributed;

res_random = res(find(res(:,1)==0), 2:8);
res_ca = res(find(res(:,1)==1), 2:8);
res_pair = res(find(res(:,1)==2), 2:8);
res_tmc = res(find(res(:,1)==4), 2:8);

f = fopen('res_taxi.txt', 'w');
for samplingi = 1 : samplingn
    sampling = samplings(samplingi);
    if (sampling == 0) 
        output_res = res_random;
        fprintf(f, 'res_random\n');
    elseif (sampling == 1)
        output_res = res_ca;
        fprintf(f, 'res_ca\n');
    elseif (sampling == 2)
        output_res = res_pair;
        fprintf(f, 'res_pair\n');
    elseif (sampling == 4)
        output_res = res_tmc;
        fprintf(f, 'res_tmc\n');
    else
        continue;
    end
    [rc,~] = size(output_res);
    for i = 1 : rc
        for j = 1 : 7
            fprintf(f, '%f ', output_res(i,j));
        end
        fprintf(f, '\n');
    end
end
fclose(f);