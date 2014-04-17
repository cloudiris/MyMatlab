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
positioncost = load('200positioncost_conn.txt');
% positioncost = ones(101,107);
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
% intvbatch = intvn;
mfrac_min = 0.1;
mfrac_gap = 0.05;
mfrac_max = 0.5;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracn = length(mfracs);

times = 1;
single_runs = mfracn * times;

samplings = [1,2];
samplingn = length(samplings);

runs = samplingn * mfracn;

f = fopen('res_taxi_intv.txt', 'w');

intvs = [1,13,26,130];
for intvi = 1 : length(intvs)
    intvbatch = intvs(intvi);
    exp_taxi_ca;

    res_random = res(find(res(:,1)==0), 2:8);
    res_ca = res(find(res(:,1)==1), 2:8);
    res_pair = res(find(res(:,1)==2), 2:8);
    res_tmc = res(find(res(:,1)==4), 2:8);

    for samplingi = 1 : samplingn
        sampling = samplings(samplingi);
        if (sampling == 0) 
            output_res = res_random;
            fprintf(f, 'res_random batch=%d\n', intvbatch);
        elseif (sampling == 1)
            output_res = res_ca;
            fprintf(f, 'res_ca batch=%d\n', intvbatch);
        elseif (sampling == 2)
            output_res = res_pair;
            fprintf(f, 'res_pair batch=%d\n', intvbatch);
        elseif (sampling == 4)
            output_res = res_tmc;
            fprintf(f, 'res_tmc batch=%d\n', intvbatch);
        else
            continue;
        end
        [rc,~] = size(output_res);
        for i = 1 : rc
            fprintf(f, '%d ', intvbatch);
            for j = 1 : 7
                fprintf(f, '%f ', output_res(i,j));
            end
            fprintf(f, '\n');
        end
    end
end
fclose(f);