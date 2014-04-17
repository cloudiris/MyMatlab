clear
clc

NDEBUG = 0;

workpath = 'd:\document\codes\matlab\wcs\';
week = 'w2';
seg = load('200seg.txt');
n = length(seg);
gt = load(['200gti_' week '.txt']);
[intvn, ~] = size(gt);
% positioncost = load('200positioncost_station.txt');
% positioncost = load('200positioncost_uniform.txt');
positioncost = load('200positioncost_conn.txt');
% positioncost = ones(101,107);
dim = size(positioncost);
% timecost = load('30Mtimecost_t_exp.txt');
% timecost = load('30Mtimecost_t_ip.txt');
% timecost = load('30Mtimecost_t_linear.txt');
% timecost = load(['30Mtimecost_' week '.txt']);
timecost = ones(9000, intvn);
recnum = load(['30Mrecnum_' week '.txt']);
tN = sum(recnum);
cabcount = load(['cabcount_' week '.txt']);
cabID = find(cabcount  > 0);
cabnum = numel(cabID);

segmapping = zeros(dim(1), dim(2));
for pt = 1 : n
    segmapping(seg(pt,1) + 1, seg(pt,2) + 1) = pt;
end

recovery = 'CS';
intvbatch = intvn;
mfrac_min = 0.05;
mfrac_gap = 0.05;
mfrac_max = 0.4;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracn = length(mfracs);

times = 1;
single_runs = mfracn * times;

samplings = [2,2.3,2.4];
samplingn = length(samplings);

runs = samplingn * mfracn;

exp_taxi_ca;

res_2 = res(find(res(:,1)==2), 2:8);
res_3 = res(find(res(:,1)==2.3), 2:8);
res_4 = res(find(res(:,1)==2.4), 2:8);

f = fopen('res_taxi_pair.txt', 'w');
for samplingi = 1 : samplingn
    sampling = samplings(samplingi);
    group = round(10*(sampling-floor(sampling)));
    if (group == 0) group = 2; end
    group
    if (group == 3) output_res = res_3;
    elseif (group == 4) output_res = res_4;
    else output_res = res_2; end
    [rc,~] = size(output_res);
    fprintf(f, 'group=%d\n', group);
    for i = 1 : rc
        for j = 1 : 7
            fprintf(f, '%f ', output_res(i,j));
        end
        fprintf(f, '\n');
    end
end
fclose(f);