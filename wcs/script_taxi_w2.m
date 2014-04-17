clear
clc

NDEBUG = 1;

workpath = 'd:\document\codes\matlab\wcs\';
seg = load('200seg.txt');
positioncost = load('200positioncost_station.txt');
timecost = load('30Mtimecost_w1.txt');
recnum = load('30Mrecnum_w1.txt');
gt = load('200gti_w1.txt');
cabcount = load('cabcount_w1.txt');
cabID = find(cabcount  > 0);
cabnum = numel(cabID);

dim = size(positioncost);
n = length(seg);
[intvn, ~] = size(gt);
tN = sum(recnum);

segmapping = zeros(dim(1), dim(2));
for pt = 1 : n
    segmapping(seg(pt,1), seg(pt,2)) = pt;
end

intvbatch = intvn;
mfrac_min = 0.05;
mfrac_gap = 0.05;
mfrac_max = 0.4;
times = 2;
mfrac = mfrac_min;

res_random = [];
res_ca = [];
res_pair = [];
res_tmc = [];

rc = 0;

for mfrac = mfrac : mfrac_gap : mfrac_max
    
    for iter = 1 : times
        rc = rc + 1;
        
        fprintf('\nRS(%.2f) ITER(%d) begin:\n', mfrac, iter);
        sampling = 0; % random
        tic;
        exp_taxi_ca;
        t = toc;
        res_random(rc, 1:7) = [afrac, acc, err, coverage, totalcost, ratiocost, t];
        fprintf('RS(%.2f) ITER(%d) finished: %f s\n', mfrac, iter, t);
        
        fprintf('\nCA(%.2f) ITER(%d) begin:\n', mfrac, iter);
        sampling = 1; % ca
        tic;
        exp_taxi_ca;
        t = toc;
        res_ca(rc, 1:7) = [afrac, acc, err, coverage, totalcost, ratiocost, t];
        fprintf('CA(%.2f) ITER(%d) finished: %f s\n', mfrac, iter, t);
        
        fprintf('\nPW(%.2f) ITER(%d) begin:\n', mfrac, iter);
        sampling = 2; % pair
        tic;
        exp_taxi_ca;
        t = toc;
        res_pair(rc, 1:7) = [afrac, acc, err, coverage, totalcost, ratiocost, t];
        fprintf('PW(%.2f) ITER(%d) finished: %f s\n', mfrac, iter, t);
        
        fprintf('\nTMC(%.2f) ITER(%d) begin:\n', mfrac, iter);
        sampling = 4; % TMC
        tic;
        exp_taxi_ca;
        t = toc;
        res_tmc(rc, 1:7) = [afrac, acc, err, coverage, totalcost, ratiocost, t];
        fprintf('TMC(%.2f) ITER(%d) finished: %f s\n', mfrac, iter, t);
    end
end

f = fopen('res_taxi.txt', 'w');
for ri = 0 : 4
    if (ri == 0) 
        res = res_random;
        fprintf(f, 'res_random\n');
    elseif (ri == 1)
        res = res_ca;
        fprintf(f, 'res_ca\n');
    elseif (ri == 2)
        res = res_pair;
        fprintf(f, 'res_pair\n');
    elseif (ri == 4)
        res = res_tmc;
        fprintf(f, 'res_tmc\n');
    else
        continue;
    end
    for i = 1 : rc
        for j = 1 : 7
            fprintf(f, '%f ', res(i,j));
        end
        fprintf(f, '\n');
    end
end
fclose(f);
    