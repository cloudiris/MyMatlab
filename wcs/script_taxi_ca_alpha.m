clear
clc

NDEBUG = 1;

workpath = 'd:\document\codes\matlab\wcs\';
week = 'w2';
seg = load('200seg.txt');
n = length(seg);
gt = load(['200gti_' week '.txt']);
[intvn, ~] = size(gt);
% positioncost = load('200positioncost_station.txt');
% positioncost = load('200positioncost_uniform.txt');
% positioncost = load('200positioncost_conn.txt');
positioncost = ones(101,107);
dim = size(positioncost);
% timecost = load('30Mtimecost_t_exp.txt');
% timecost = load('30Mtimecost_t_ip.txt');
timecost = load('30Mtimecost_t_linear.txt');
% timecost = load(['30Mtimecost_' week '.txt']);
% timecost = ones(9000, intvn);
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
mfrac_min = 0.1;
mfrac_gap = 0.05;
mfrac_max = 0.7;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracn = length(mfracs);

times = 1;
single_runs = mfracn * times;

% samplings = [1,1.1,1.2,1.3,1.4,1.5];
samplings = [1,1.1,1.2,1.3];
samplingn = length(samplings);

runs = samplingn * mfracn;

exp_taxi_ca;

res_1 = res(find(res(:,1)==1), 2:8);
res_dot1 = res(find(res(:,1)==1.1), 2:8);
res_dot5 = res(find(res(:,1)==1.2), 2:8);
res_2 = res(find(res(:,1)==1.3), 2:8);
res_5 = res(find(res(:,1)==1.4), 2:8);
res_100 = res(find(res(:,1)==1.5), 2:8);

f = fopen('res_taxi_ca.txt', 'w');
for samplingi = 1 : samplingn
    sampling = samplings(samplingi);
    if (sampling == 1) 
        output_res = res_1;
        alpha = 1;
    elseif (sampling == 1.1)
        output_res = res_dot1;
        alpha = 0.1;
    elseif (sampling == 1.2)
        output_res = res_dot5;
        alpha = 0.5;
    elseif (sampling == 1.3)
        output_res = res_2;
        alpha = 2;
    elseif (sampling == 1.4)
        output_res = res_5;
        alpha = 5;
    elseif (sampling == 1.5)
        output_res = res_100;
        alpha = 100;
    else
        continue;
    end
    [rc,~] = size(output_res);
    fprintf(f, 'alpha=%.1f\n', alpha);
    for i = 1 : rc
        fprintf(f, '%f ', alpha);
        for j = 1 : 7
            fprintf(f, '%f ', output_res(i,j));
        end
        fprintf(f, '\n');
    end
end
fclose(f);