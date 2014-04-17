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
% positioncost = load('test_positioncost.txt');
% positioncost = load('200positioncost_uniform.txt');
positioncost = ones(101,107);
dim = size(positioncost);
% timecost = ones(1, intvn);
timecost = load(['30Mtimecost_' week '.txt']);
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

samplings = [0,1,2,4];
samplingn = length(samplings);

runs = samplingn * mfracn;

f = fopen('res_taxi_timenoise.txt', 'w');

noiselevel = [0.01,0.1,1,10,100];
timecost_bk = timecost;
for noisei = 1 : 5
    tnoise = noise_n(intvn, sqrt(noiselevel(noisei)*norm(timecost_bk)^2))';
    timecost = abs(timecost_bk + tnoise);
    anoise = 20*log10(norm(timecost-timecost_bk)/norm(timecost_bk));
    exp_taxi_ca;

    for samplingi = 1 : 4
        sampling = samplings(samplingi);
        res_random = res(find(res(:,1)==0), 2:8);
        res_ca = res(find(res(:,1)==1), 2:8);
        res_pair = res(find(res(:,1)==2), 2:8);
        res_tmc = res(find(res(:,1)==4), 2:8);

        if (sampling == 0) 
            output_res = res_random;
            fprintf(f, 'res_random noise=%f\n', anoise);
        elseif (sampling == 1)
            output_res = res_ca;
            fprintf(f, 'res_ca noise=%f\n', anoise);
        elseif (sampling == 2)
            output_res = res_pair;
            fprintf(f, 'res_pair noise=%f\n', anoise);
        elseif (sampling == 4)
            output_res = res_tmc;
            fprintf(f, 'res_tmc noise=%f\n', anoise);
        else
            continue;
        end
        [rc,~] = size(output_res);
        for i = 1 : rc
            fprintf(f, '%f ', anoise);
            for j = 1 : 7
                fprintf(f, '%f ', output_res(i,j));
            end
            fprintf(f, '\n');
        end
    end
end
fclose(f);